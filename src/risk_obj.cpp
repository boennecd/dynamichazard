#include <iostream>
#include <Rcpp.h>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]
using namespace std;
using namespace Rcpp;



#if defined(NDEBUG)
#undef NDEBUG
#endif


// use vecmin from other cpp function
extern int vecmin(IntegerVector x);

//' Function to find the risk set in c++
// [[Rcpp::export]]
List get_risk_obj_rcpp(const NumericVector &start, const NumericVector &stop,
                       const LogicalVector &event,
                       const double &by, // -1 if missing
                       const IntegerVector &start_order, // index of stop sorted in decreasing order
                       const double &max_T,
                       const IntegerVector &order_by_id_and_rev_start,
                       const IntegerVector &id
){
  // See this page for rcpp sugar functions http://dirk.eddelbuettel.com/code/rcpp/html/unique_8h.html
  NumericVector stop_events = stop[event];
  stop_events = stop_events[stop_events <= max_T];

  const double min_start = min(start);

  double I_stop_time, I_start_time;
  const unsigned int n = start.size();
  unsigned int d, i, j, d_;
  NumericVector event_times;
  vector<int> indicies; // see http://stackoverflow.com/questions/13324431/c-vectors-insert-push-back-difference

  // Find event times
  if(by < 0.0){
    event_times = sort_unique(stop_events);
    d = event_times.size();
  }
  else{
    d = ceil((max_T - min_start) / by);
    event_times = NumericVector(d);
    for(i = 0; i < d; i++){
      event_times[i] = min_start + (i + 1) * by;
    }
  }

  // Find new rounded stop times
  NumericVector stop_new = clone(stop);
  I_start_time = min_start;
  for(d_ = 0; d_ < d + 1; d_++){
    if(d_ == d)
      I_stop_time = max(stop) + 1.0;
    else
      I_stop_time = event_times[d_];

    for(i = 0; i < n; i++){
      j = start_order[i];

      if(start[j] >= I_stop_time)
        break;

      if(I_start_time < stop[j] && stop[j] <= I_stop_time){
        stop_new[j] = I_stop_time;
      }
    }

    I_start_time = I_stop_time;
  }

  // Make flag for whether an event is inside a bin
  LogicalVector is_within_bin_and_event(n);
  LogicalVector new_events_flags = clone(event);

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

      if(start[i] >= I_stop_time)
        break;

      if(start[i] <= I_start_time && I_start_time < stop[i])
        indicies.push_back(i + 1); // R use one indexing!
      else if(start[i] >= I_start_time && stop[i] <= I_stop_time)
        is_within_bin_and_event[i] = event[i];
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

  unsigned int n_events_inside = 0;
  int this_id, this_bin;
  unsigned int k, l;
  double bin_start, bin_stop;
  double* tmp_pointer;

  j = 0;
  while(j < n){
    i = order_by_id_and_rev_start[j];

    // Skipe those that are not events insides bins
    if(!is_within_bin_and_event[i]){
      j++;
      continue;
    }

    n_events_inside++;

    // Find the bin (see http://stackoverflow.com/a/15724226/5861244)
    tmp_pointer =  lower_bound(dum_for_find_interval.begin(), dum_for_find_interval.end(), stop[i]);
    this_bin = tmp_pointer - dum_for_find_interval.begin() - 1;
    this_id = id[i];

    bin_start = dum_start[this_bin];
    bin_stop = dum_stop[this_bin];

    if(j == n - 1){
      j++;
      continue;
    }

    for(k = j + 1; k < n; k++){
      l = order_by_id_and_rev_start[k];

      if(id[l] != this_id ||
         stop[l] <= bin_start){
        j = k;
        break;
      }
      else if(start[l] <= bin_start){
        new_events_flags[l] = true;
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
                      Named("event_times") = event_times,
                      Named("I_len") = I_len,
                      Named("d") = d,
                      Named("stop_new") = stop_new,
                      Named("new_events_flags") = new_events_flags));

}
