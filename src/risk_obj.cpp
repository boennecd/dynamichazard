#include "dynamichazard.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List get_risk_obj_rcpp(const NumericVector &start, const NumericVector &stop,
                       const LogicalVector &event,
                       const double &by, // -1 if missing
                       const IntegerVector &start_order, // index of stop sorted in decreasing order
                       const double &max_T,
                       const IntegerVector &order_by_id_and_rev_start,
                       const IntegerVector &id,
                       const bool &is_for_discrete_model = true
){
  // See this page for rcpp sugar functions http://dirk.eddelbuettel.com/code/rcpp/html/unique_8h.html
  NumericVector stop_events = stop[event];
  stop_events = stop_events[stop_events <= max_T];

  const double delta_t_eps = 1.0e-14;

  const double min_start = min(start);

  double I_stop_time, I_start_time; // current interval start and stop time
  const unsigned int n = start.size();
  unsigned int d, i, j, d_;
  NumericVector event_times; // vector of event times
  vector<int> indicies; // see http://stackoverflow.com/questions/13324431/c-vectors-insert-push-back-difference

  // Find event times
  if(by < 0.0){ // set event times to all the unique stop times
    event_times = sort_unique(stop_events);
    d = event_times.size();
  }
  else{ // go from min to max with by as increaments
    d = ceil((max_T - min_start) / by);
    event_times = NumericVector(d);
    for(i = 0; i < d; i++){
      event_times[i] = min_start + (i + 1) * by;
    }
  }

  // Make flag for whether an observation is inside a bin
  //  E.g. we have a row at time [.25, .75) with an event
  //  and a bins with times [0, 1).
  LogicalVector is_within_bin(n);
  //LogicalVector new_events_flags = clone(event); // TODO: Delete

  // Find risket sets
  List risk_sets(d);
  I_stop_time = event_times(d - 1);

  // Set bin start time
  for(d_ = d - 1; ; d_--){
    if(d_ == 0)
      I_start_time = min_start;
    else
      I_start_time = event_times(d_ - 1);

    indicies = {};
    for(j = 0; j < n; j++){
      i = start_order[j];

      if(start[i] >= I_stop_time) // no need to search further
        break;

      if(is_for_discrete_model){
        if(start[i] <= I_start_time && I_start_time < stop[i]){
          // starts before and ends after. Thus in this bin
          indicies.push_back(i + 1); // R use one indexing!
        }
        else if(start[i] >= I_start_time && stop[i] <= I_stop_time){
          // start after and ends before. Thus inside the bin
          is_within_bin[i] = true;
        }
      } else {
        if(I_start_time < stop[i]){
          // we all ready know it start before the interval stops and now we
          // know it stop after the interval starts so it is in
          /*if(std::min(stop[i], I_stop_time) - std::max(start[i], I_start_time) < delta_t_eps){ //TODO: Delete
            // there is numerical imprecission with float comparisson. Thus the check is added
            // TODO: is there something better?
            continue;
          }*/

          indicies.push_back(i + 1); // R use one indexing!
        }
      }
    }

    risk_sets[d_] = IntegerVector(indicies.begin(), indicies.end());

    I_stop_time = I_start_time;

    if(d_ == 0)
      break;
  }

  // Hanlde events inside bins
  NumericVector dum_for_find_interval(d + 2);
  NumericVector dum_start(d + 1);
  NumericVector dum_stop(d + 1);

  dum_start[0] = dum_for_find_interval[0] = min_start;

  for(i = 1; i < d + 1; i++){
    dum_start[i] = dum_stop[i - 1] = dum_for_find_interval[i] = event_times[i - 1];
  }

  dum_for_find_interval[d + 1] = dum_stop[d] = max(max(stop) + 1.0, dum_stop[d]);

  int this_id, this_bin;
  unsigned int k, l;
  double bin_start, bin_stop;
  double *tmp_pointer;

  // this vector is used to indicate if an observation is an event in a given
  // bin. The default is -1 which indicates that the observation is not an
  // event in any bin. Note that this is zero indexed
  IntegerVector is_event_in(n, -1);

  j = 0;
  while(j < n){
    i = order_by_id_and_rev_start[j];

    if(!event[i]){
      // The default indicating no event is correct. Thus, we continue
      j++;
      continue;
    }

    // Find the bin (see http://stackoverflow.com/a/15724226/5861244)
    tmp_pointer =  lower_bound(dum_for_find_interval.begin(), dum_for_find_interval.end(), stop[i]);
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

    this_bin = tmp_pointer - dum_for_find_interval.begin() - 1;
    this_id = id[i];

    bin_start = dum_start[this_bin];
    bin_stop = dum_stop[this_bin];

    if(j == n - 1){ // we reached the end. No further work needed
      j++;
      break;
    }

    for(k = j + 1; k < n; k++){
      l = order_by_id_and_rev_start[k];

      if(id[l] != this_id ||
         stop[l] <= bin_start){
        // Either there is a gab for this id or we changed to another indvidual
        // we break in either case
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
  NumericVector temp = event_times;
  temp.push_front(min_start);
  NumericVector I_len = diff(temp);

  return(List::create(Named("risk_sets") = risk_sets,
                      Named("min_start") = min_start,
                      Named("event_times") = event_times,
                      Named("I_len") = I_len,
                      Named("d") = d,
                      Named("is_event_in") = is_event_in,
                      Named("is_for_discrete_model") = is_for_discrete_model));
}
