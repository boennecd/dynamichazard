template<typename iterator>
  class work_block {
    public:
      iterator start;
    iterator end;
    const unsigned long block_size;

    work_block(iterator start, iterator end, const unsigned long block_size):
      start(start), end(end), block_size(block_size) {}
  };

template<typename iterator>
  std::vector<work_block<iterator>> get_work_blocks(
    iterator begin, iterator end, const unsigned long block_size){
    unsigned long const length = std::distance(begin, end);
    unsigned long const num_blocks= (length + block_size - 1) / block_size;

    std::vector<work_block<iterator>> ans;
    ans.reserve(num_blocks);

    iterator block_start = begin;
    for(unsigned long i = 0; i < num_blocks - 1; ++i){
      iterator block_end = block_start;
      std::advance(block_end, block_size);

      ans.emplace_back(block_start, block_end, block_size);
      block_start = block_end;
    }

    ans.emplace_back(block_start, end, std::distance(block_start, end));

    return ans;
  }
