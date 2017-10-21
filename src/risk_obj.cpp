#include "arma_n_rcpp.h"

#include <map>

// [[Rcpp::export]]
Rcpp::List get_risk_obj_rcpp(
    const Rcpp::NumericVector &start, const Rcpp::NumericVector &stop,
    const Rcpp::LogicalVector &event,
    const double &by, // -1 if missing
    const Rcpp::IntegerVector &start_order, // index of stop sorted in decreasing order
    const double &max_T,
    const Rcpp::IntegerVector &order_by_id_and_rev_start,
    const Rcpp::IntegerVector &id,
    const double min_start,
    const Rcpp::NumericVector &event_times_in,
    const bool &is_for_discrete_model = true){
  double I_stop_time, I_start_time; // current interval start and stop time
  const unsigned int n = start_order.size(); // Number of observations
  unsigned int i, j, d_;
  Rcpp::NumericVector event_times = Rcpp::clone(event_times_in);
  unsigned int d = event_times.size();
  std::vector<int> indicies;

  if(start_order.size() != order_by_id_and_rev_start.size())
    Rcpp::stop("Size of start_order and order_by_id_and_rev_start does match");

  // Make flag for whether an observation is inside a bin
  //  E.g. we have a row at time [.25, .75) and a bin
  //  with times [0, 1) then this observation is inside a
  //  bin.
  // The vector will be false if is_for_discrete_model is
  // true
  // The length has to match start, stop, id but not
  // order_by_id_and_rev_start and start_order
  Rcpp::LogicalVector is_within_bin(start.size());

  // List of indicies in the risk sets in each of the d bins
  Rcpp::List risk_sets(d);

  I_stop_time = event_times(d - 1);
  for(d_ = d - 1; ; d_--){ // Loop over intervals
    if(d_ == 0)
      I_start_time = min_start;
    else
      I_start_time = event_times(d_ - 1);

    indicies = {};
    for(j = 0; j < n; j++){ // Loop over observations
      i = start_order[j];

      if(start[i] >= I_stop_time) // no need to search further due to ordering
        break;

      if(is_for_discrete_model){
        if(start[i] <= I_start_time && I_start_time < stop[i]){
          // Starts before and ends after start
          //   E.g. the bin is [1, 2) then these observation these are in:
          //    * [0.5, 1.5) and [1, 1.5)
          //   These are not in:
          //    * [1.1, 1.5) and [0.5, 1)
          indicies.push_back(i + 1); // R use one indexing!
        }
        else if(start[i] > I_start_time && stop[i] <= I_stop_time){
          // Start after start time and end before stop time
          //   E.g. the bin is [1, 2) then this observation is inside:
          //    *  [1.1, 1.5)
          //   These ones are not
          //    *  [0.5, 2) and [1.5, 2.5)
          is_within_bin[i] = true;
        }
      } else {
        if(I_start_time < stop[i]){
          // Start before stop time (we know this due to the former check with
          // the break statement) and ends after the start time
          indicies.push_back(i + 1); // R use one indexing!
        }
      }
    }

    risk_sets[d_] = std::move(
      Rcpp::IntegerVector(indicies.begin(), indicies.end()));

    I_stop_time = I_start_time;

    if(d_ == 0)
      break;
  }

  // Set the is_event_in flag and if is_for_discrete_model is true then remove
  // some item from the risk set that stops inside a bin with right-censoring
  // An example of the later is if we have an inividual with times [0, 1) and
  // [1, 1.5) with no event at the end the bin are [0, 1) and [1, 2). We do not
  // know if the indvidual survives [1.5, 2) so we remove his observation from
  // bin 2
  Rcpp::NumericVector dum_for_find_interval(d + 2);
  Rcpp::NumericVector dum_start(d + 1);
  Rcpp::NumericVector dum_stop(d + 1);

  dum_start[0] = dum_for_find_interval[0] = min_start;

  for(i = 1; i < d + 1; i++){
    dum_start[i] = dum_stop[i - 1] = dum_for_find_interval[i] = event_times[i - 1];
  }

  dum_for_find_interval[d + 1] = dum_stop[d] =
    std::max(max(stop) + 1.0, dum_stop[d]);

  int old_id, this_id, this_bin;
  unsigned int k, l;
  double bin_start, bin_stop;
  double *tmp_pointer;

  // This vector is used to indicate if an observation is an event in a given
  // bin. The default is -1 which indicates that the observation is not an
  // event in any bins. Note that this is zero indexed as it will only be used
  // in cpp
  // The length has to match start, stop, id but not
  // order_by_id_and_rev_start and start_order
  Rcpp::IntegerVector is_event_in(start.size(), -1);

  j = 0;
  this_id = std::numeric_limits<int>::min(); // NB: Assuming this not one of
                                             // ids!
  while(j < n){ // Loop over rows in id and decreasing start time order
    old_id = this_id;
    i = order_by_id_and_rev_start[j];

    this_id = id[i];

    if(!event[i]){
      // The default indicating no event is correct
      if(is_for_discrete_model && this_id != old_id){
        // Remove observation that are right censored ending inside a bin from
        // the risk set
        //   E.g we have an inividual with times [0, 1) and [1, 1.5) with no
        //   event at the end the bin are [0, 1) and [1, 2). We do not know if
        //   the indvidual survives [1.5, 2) so we remove his observation from
        //   bin 2
        //   Another example is an inividual with observations [0, 1.5),
        //   [1.5, 1.6) and [1.6, 1.8). Here we need to remove [0, 1.5) from
        //   risk set for bin [1, 2)
        // This can be done by:
        //  1. Find the bin number
        //  2. Check if the stop times is strictly before the bin stop time
        //  3. Remove the observation for this id that is inside the bin found
        //     in 1. from the corresponding risk set

        // Find the bin that the last observation stops in
        tmp_pointer =  std::lower_bound(
          dum_for_find_interval.begin(), dum_for_find_interval.end(), stop[i]);
        int bin_number = std::distance(
          dum_for_find_interval.begin(), tmp_pointer) - 1;

        bin_stop = dum_stop[bin_number];

        if(stop[i] == bin_stop){
          // No need to do any thing further. E.g. we have an observation
          // [1, 2) with bin [1, 2)
          j++;
          continue;
        }

        bin_start = dum_start[bin_number];

        // All observation from this point on will have stop[i] < bin_stop
        // Thus, we need to fint the first with:
        //  * bin_start < stop[i] and start[i] <= bin_start and remove it from
        //    the risk set or ...
        //  * stop[i] <= bin_start in which case we need to break or
        //  * The id changed in which case we can continue

        if(bin_number < event_times.size()){
          // The last observation ends before the last bin so further action is
          // required
          for(k = j; k < n; k++){
            l = order_by_id_and_rev_start[k];

            if(this_id != id[l]){
              // Id changed so we continue
              j = k;
              break;
            }

            if(bin_start < stop[l] && start[l] <= bin_start){
              // Needs to be removed. This could properly be written smarter
              // TODO: do this smarter
              std::vector<int> r_set = Rcpp::as<std::vector<int>>(risk_sets[bin_number]);
              r_set.erase(std::remove(r_set.begin(), r_set.end(), l + 1), // we entered the value as one indexed
                          r_set.end());
              risk_sets[bin_number] = Rcpp::IntegerVector(r_set.begin(), r_set.end());

              j = k;
              break;

            } else if(stop[l] <= bin_start){
              j = k;
              break;
            }
          }
        }

        continue;

      } else{

        j++;
        continue;
      }
    }

    // Find the bin that the last observation stops in
    tmp_pointer =  std::lower_bound(dum_for_find_interval.begin(), dum_for_find_interval.end(), stop[i]);
    int bin_number = std::distance(dum_for_find_interval.begin(), tmp_pointer) - 1;

    if(bin_number == event_times.size()){
      // The failure is after the last event time
      // Thus, we do not include the failure
      j++;
      continue;
    }

    if(!is_within_bin[i]){
      // This an event and it do cross two bins. Thus, we need to set the bin
      // number in the is_event_in
      is_event_in[i] = bin_number;
      j++;
      continue;
    }

    // Here we need to handle the case where the observation does not cross two
    // bins
    //   E.g. we have two observations [0, 1.5) and [1.6, 1,8) the later with
    //   an event with bins [0, 1) and [1, 2)
    //   Then the [0, 1.5) should be marked as an event in [1, 2)
    this_bin = tmp_pointer - dum_for_find_interval.begin() - 1;
    bin_start = dum_start[this_bin];

    if(j == n - 1){ // we reached the end. No further work needed
      j++;
      break;
    }

    for(k = j + 1; k < n; k++){
      l = order_by_id_and_rev_start[k];

      if(id[l] != this_id ||
         stop[l] <= bin_start){
        // Either there is a gab for this id. E.g. we have observations
        // [0, 0.8) and [1.2, 1.5) with bins [0, 1) and [1, 2) or the id have
        // changed. In both cases we break
        j = k;
        break;
      }
      else if(start[l] <= bin_start){
        // We found the previous row. We mark it as an event and continue
        is_event_in[l] = bin_number;
        j = k;
        break;
      }

      j = k;
    }
  }

  // compute length of each interval
  event_times.push_front(min_start);
  Rcpp::NumericVector I_len = diff(event_times);

  return(Rcpp::List::create(Rcpp::Named("risk_sets") = risk_sets,
                            Rcpp::Named("min_start") = min_start,
                            Rcpp::Named("event_times") = event_times,
                            Rcpp::Named("I_len") = I_len,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("is_event_in") = is_event_in,
                            Rcpp::Named("is_for_discrete_model") = is_for_discrete_model));
}
