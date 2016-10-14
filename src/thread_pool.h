#include <deque>
#include <future>
#include <memory>
#include <functional>
#include <iostream>
#include <iostream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <type_traits>

class thread_pool
{
public:
	thread_pool(unsigned int n_jobs);

	template<typename FunctionType>
	std::future<typename std::result_of<FunctionType()>::type> submit(FunctionType f);

	void run_pending_task();
};
