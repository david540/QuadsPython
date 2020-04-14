class DisjointSetWithSign {

public:

  // constructor knowing the number of elements

  explicit DisjointSetWithSign(size_t num)

    : m_ids()

    , m_signs()

  {

    m_ids.resize(num);

    for (index_t i = 0; i < index_t(num); i++)

      m_ids[i] = i;

    m_signs.resize(num, true);

  }

  DisjointSetWithSign(DisjointSetWithSign const &)=delete;

  DisjointSetWithSign &operator=(DisjointSetWithSign const &)=delete;

  // merge two component

  void merge(int a, int b, bool sameSign){

    geo_assert(a>=0 && a<int(m_ids.size()));

    geo_assert(b>=0 && b<int(m_ids.size()));

    connectToRoot(index_t(a));

    connectToRoot(index_t(b));

    m_ids[m_ids[index_t(a)]] = m_ids[index_t(b)];

    m_signs[m_ids[index_t(a)]]=((m_signs[index_t(a)] == m_signs[index_t(b)]) == sameSign);

  }

​

  //

  // access to data

  //



  // return the number of sets

  int getNumSets() const {

    // connect the root

    for (index_t i = 0; i < m_ids.size(); i++) connectToRoot(i);

    // count the different roots

    int num = 0;

    for (index_t i = 0; i < m_ids.size(); i++) {

      if (i == m_ids[i]) ++num;

    }

    return num;

  }

  // return the number of set and fill setId

  template <class T>

  int getSetsId(vector<T> &idToSetId) const {

    // connect the root

    for (index_t i = 0; i < m_ids.size(); i++) connectToRoot(i);

    // prepare the correspondance root Id => setId

    idToSetId.resize(m_ids.size());

    index_t numSets = 0;

    for (index_t i = 0; i < m_ids.size(); i++) {

      if (i == m_ids[i]) {

        idToSetId[i] = numSets;

        numSets++;

      }

    }

    // fill the other set Ids

    for (index_t i = 0; i < m_ids.size(); i++)

      idToSetId[T(i)] = idToSetId[T(m_ids[i])];

    return int(numSets);

  }

​

  // return the number of values

  size_t size() const {

    return m_ids.size();

  }

  // return the parent indices

  vector<index_t> const &data() const {

    return m_ids;

  }

  /* return the signs

​

     \note must be called/used after calling getSetsId and before doing any merge */

  vector<bool> const &signs() const {

    return m_signs;

  }

protected:

  void connectToRoot(index_t i) const {

    // find the root

    index_t rootId=m_ids[i];

    while (rootId!=m_ids[rootId]) rootId=m_ids[rootId];

    // connect all the path to root

    while (i!=rootId) {

      auto newI=m_ids[i];

      m_ids[i] = rootId;

      m_signs[i] = m_signs[i] == m_signs[rootId];

      i=newI;

    }

  }

  mutable vector<index_t> m_ids;

  mutable vector<bool> m_signs;

};
