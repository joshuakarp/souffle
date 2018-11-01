#pragma once

#include <iostream>
#include <atomic>
#include <cstring>
#include <array>
#include <list>
#include "ParallelUtils.h"

using std::size_t;
namespace souffle {

/**
 * A PiggyList that allows insertAt functionality.
 * This means we can't append, as we don't know the next available element.
 * insertAt is dangerous. You must be careful not to call it for the same index twice!
 */
template <class T>
class RandomInsertPiggyList {
    const size_t BLOCKBITS = 16ul;
    const size_t INITIALBLOCKSIZE = (1ul << BLOCKBITS);

    // number of elements currently stored within
    std::atomic<size_t> numElements;

    // 2^64 - 1 elements can be stored (default initialised to nullptrs)
    static constexpr size_t maxContainers = 64;
    std::array<std::atomic<T*>, maxContainers> blockLookupTable = {};
    
    // for parallel node insertions
    mutable SpinLock slock;

    /**
     * Free the arrays allocated within the linked list nodes
     */
    void freeList() {
        // delete all - deleting a nullptr is a no-op
        for (size_t i = 0; i < maxContainers; ++i) {
            delete[] blockLookupTable[i].load();
            // reset the container within to be empty.
            blockLookupTable[i].store(nullptr);
        }
    }

    public:
    RandomInsertPiggyList() : numElements(0) { }
    // an instance where the initial size is not 65k, and instead is user settable (to a power of initialbitsize)
    RandomInsertPiggyList(size_t initialbitsize) : BLOCKBITS(initialbitsize), numElements(0) {}

    /** copy constructor */
    RandomInsertPiggyList(const RandomInsertPiggyList& other) : BLOCKBITS(other.BLOCKBITS) {
        this->numElements.store(other.numElements.load());

        // copy blocks from the old lookup table to this one
        for (size_t i = 0; i < maxContainers; ++i) {
            if (other.blockLookupTable[i].load() != nullptr) {
                // calculate the size of that block
                const size_t blockSize = INITIALBLOCKSIZE << i;

                // allocate that in the new container
                this->blockLookupTable[i].store(new T[blockSize]);

                // then copy the stuff over
                std::memcpy(this->blockLookupTable[i].load(), other.blockLookupTable[i].load(), blockSize * sizeof(T));
            }
        }
    }

    // TODO: implement these ctors (...maybe not the move ones..?)
    RandomInsertPiggyList(RandomInsertPiggyList&& other) = delete;
    /** move constructor (fastest move ctor is a copy, due to atomics) */
    //RandomInsertPiggyList(RandomInsertPiggyList&& other) : BLOCKBITS(other.BLOCKBITS) {
        //this->numElements.store(other.numElements.load());

        //for (size_t i = 0; i < maxContainers; ++i) {
            //this->blockLookupTable[i].store(other.blockLookupTable[i].load());
        //}
    //}

    // copy assign ctor (can't actually be legit done, as we contain const data)
    RandomInsertPiggyList& operator=(RandomInsertPiggyList& other) = delete;
    //RandomInsertPiggyList& operator=(RandomInsertPiggyList& other) {
    //    if (&other == this) return *this;

    //    freeList();

    //    this->numElements.store(other.numElements.load());
    //    for (size_t i = 0; i < maxContainers; ++i) {
    //        this->blockLookupTable[i].store(other.blockLookupTable[i].load());
    //    }
    //    return *this;
    //}

    // move assign ctor
    RandomInsertPiggyList& operator=(RandomInsertPiggyList&& other) = delete;

    ~RandomInsertPiggyList() {
        freeList();
    }

    inline size_t size() const { return numElements.load(); }

    inline T* getBlock(size_t blockNum) const { return blockLookupTable[blockNum]; }

    inline T& get(size_t index) const {
        size_t nindex = index + INITIALBLOCKSIZE;
        size_t blockNum = (63 - __builtin_clzll(nindex));
        size_t blockInd = (nindex) & ((1 << blockNum) - 1);
        return this->getBlock(blockNum-BLOCKBITS)[blockInd];
    }

    void insertAt(size_t index, T value) {
        // starting with an initial blocksize requires some shifting to transform into a nice powers of two series
        size_t blockNum = (63 - __builtin_clzll(index + INITIALBLOCKSIZE)) - BLOCKBITS;

        // allocate the block if not allocated
        if (blockLookupTable[blockNum].load(std::memory_order_relaxed) == nullptr) {
            slock.lock();
            if (blockLookupTable[blockNum].load(std::memory_order_relaxed) == nullptr) {
                blockLookupTable[blockNum].store(new T[INITIALBLOCKSIZE << blockNum]);
            }
            slock.unlock();
        }

        this->get(index) = value;
        // we ALWAYS increment size, even if there was something there before (its impossible to tell!)
        // the onus is up to the user to not call this for an index twice
        ++numElements;
    }

    void clear() {
        freeList();
        numElements.store(0);
    }
};

template <class T>
class PiggyList {
    const size_t BLOCKBITS = 16ul;
    const size_t BLOCKSIZE = (1ul << BLOCKBITS);

    // number of inserted 
    std::atomic<size_t> num_containers;
    size_t allocsize = BLOCKSIZE;
    std::atomic<size_t> container_size;
    std::atomic<size_t> m_size;

    // > 2^64 elements can be stored (default initialise to nullptrs)
    static constexpr size_t max_conts = 64;
    std::array<T*, max_conts> blockLookupTable = {};

    // for parallel node insertions
    mutable SpinLock sl;

    /**
     * Free the arrays allocated within the linked list nodes
     */
    void freeList() {
        // we don't know which ones are taken up!
        for (size_t i = 0; i < num_containers; ++i) {
            delete[] blockLookupTable[i];
        }
    }

public:
    PiggyList() : num_containers(0), container_size(0), m_size(0)  { }
    PiggyList(size_t initialbitsize) : BLOCKBITS(initialbitsize), num_containers(0), container_size(0), m_size(0) {}

    /** copy constructor */
    PiggyList(const PiggyList& other) : BLOCKBITS(other.BLOCKBITS) {

        num_containers.store(other.num_containers.load());
        container_size.store(other.container_size.load());
        m_size.store(other.m_size.load());
        // copy each chunk from other into this
        // the size of the next container to allocate
        size_t cSize = BLOCKSIZE;
        for (size_t i = 0; i < other.num_containers; ++i) {
            this->blockLookupTable[i] = new T[cSize];
            std::memcpy(this->blockLookupTable[i], other.blockLookupTable[i], cSize*sizeof(T));
            cSize <<= 1;
        }
        // if this isn't the case, uhh
        assert((cSize >> 1) == container_size.load());
    }

    /** move constructor */
    PiggyList(PiggyList&& other) = delete;
    /** copy assign ctor **/
    PiggyList& operator=(const PiggyList& other) = delete;

    ~PiggyList() {
        freeList();
    }

    /**
     * Well, returns the number of nodes exist within the list + number of nodes queued to be inserted
     *  The reason for this, is that there may be many nodes queued up
     *  that haven't had time to had containers created and updated
     * @return the number of nodes exist within the list + number of nodes queued to be inserted
     */
    inline size_t size() const {
        return m_size.load();
    };
    

    inline T* getBlock(size_t blocknum) const {
        return this->blockLookupTable[blocknum];
    }

    size_t append(T element) {
        size_t new_index = m_size.fetch_add(1, std::memory_order_acquire);

        // will this not fit?
        if (container_size < new_index + 1) {
            sl.lock();
            // check and add as many containers as required
            while (container_size < new_index + 1) {
                blockLookupTable[num_containers] = new T[allocsize];
                num_containers += 1;
                container_size += allocsize;
                // double the number elements that will be allocated next time
                allocsize <<= 1;
            }
            sl.unlock();
        }
        
        this->get(new_index) = element;
        return new_index;
    }

    size_t createNode() {
        size_t new_index = m_size.fetch_add(1, std::memory_order_acquire);

        // will this not fit?
        if (container_size < new_index + 1) {
            sl.lock();
            // check and add as many containers as required
            while (container_size < new_index + 1) {
                blockLookupTable[num_containers] = new T[allocsize];
                num_containers += 1;
                container_size += allocsize;
                // double the number elements that will be allocated next time
                allocsize <<= 1;
            }
            sl.unlock();
        }
        
        return new_index;
    }

    /**
     * Retrieve a reference to the stored value at index
     * @param index position to search
     * @return the value at index
     */
    // XXX: in another commit i made this const size_t&  - is this really that beneficial?
    inline T& get(size_t index) const {
        // supa fast 2^16 size first block
        size_t nindex = index + BLOCKSIZE;
        size_t blockNum = (63 - __builtin_clzll(nindex));
        size_t blockInd = (nindex) & ((1 << blockNum) - 1);
        return this->getBlock(blockNum-BLOCKBITS)[blockInd];
    }

    /**
     * Clear all elements from the PiggyList
     */
    void clear() {
        freeList();
        m_size = 0;
        num_containers = 0;
        
        allocsize = BLOCKSIZE;
        container_size = 0;
    }

    class iterator : std::iterator<std::forward_iterator_tag, T> {
        size_t cIndex = 0;
        PiggyList* bl;

    public:
        // default ctor, to silence
        iterator(){};

        /* begin iterator for iterating over all elements */
        iterator(PiggyList* bl) : bl(bl){};
        /* ender iterator for marking the end of the iteration */
        iterator(PiggyList* bl, size_t beginInd) : cIndex(beginInd), bl(bl){};

        T operator*() {
            return bl->get(cIndex);
        };
        const T operator*() const {
            return bl->get(cIndex);
        };

        iterator& operator++(int) {
            ++cIndex;
            return *this;
        };

        iterator operator++() {
            iterator ret(*this);
            ++cIndex;
            return ret;
        };

        friend bool operator==(const iterator& x, const iterator& y) {
            return x.cIndex == y.cIndex && x.bl == y.bl;
        };

        friend bool operator!=(const iterator& x, const iterator& y) {
            return !(x == y);
        };
    };

    iterator begin() {
        return iterator(this);
    };
    iterator end() {
        return iterator(this, size());
    };
};

}
