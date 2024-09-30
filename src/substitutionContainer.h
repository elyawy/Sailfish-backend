#ifndef _SUBSTITUTION_CONTAINER_H_
#define _SUBSTITUTION_CONTAINER_H_

#include <memory>
#include <vector>
#include <string>
#include <tuple>
#include <array>

#include <map>

#include "../libs/Phylolib/includes/sequenceContainer.h"

#define CHUNK_SIZE 2048

// #define MAX_COMPRESS_SIZE (size_t)CHUNK_SIZE/3


typedef std::map<size_t, std::vector<ALPHACHAR>> compressedChunk;

class substitutionContainer
{
private:
    std::shared_ptr<std::vector<compressedChunk>> _substitutionChunks;
    std::shared_ptr<std::vector<size_t>> _rowIndexToId;
    std::shared_ptr<std::vector<std::string>> _rowIndexToName;
    const alphabet* _alphabet;

public:
    substitutionContainer(const alphabet *alph) {
        _substitutionChunks = std::make_shared<std::vector<compressedChunk>>();

        _rowIndexToId = std::make_shared<std::vector<size_t>>();
        _rowIndexToName = std::make_shared<std::vector<std::string>>();
        _alphabet = alph->clone();

    }

    void addChunk(const sequenceContainer &seqContainer) {
        
        compressedChunk chunk;

        for (sequenceContainer::constTaxaIterator taxa = seqContainer.constTaxaBegin(); 
            taxa != seqContainer.constTaxaEnd(); ++taxa) {
            size_t id = taxa->id();
            std::string seqName = taxa->name();
            std::string currentSeqStr = taxa->toString();
            size_t currentSeqLength = taxa->seqLen();

            size_t originalSequenceLength = currentSeqStr.size();
            std::vector<ALPHACHAR> encodedBitSequence;
            // handle remainder sequence
            if (currentSeqLength != CHUNK_SIZE) {


            }

            if (_alphabet->size() == 4) encodedBitSequence.resize(CHUNK_SIZE / 4);
            if (_alphabet->size() == 20) encodedBitSequence.resize(5 * CHUNK_SIZE / 8);

            // const char* currentSeqPtr = currentSeqStr.c_str();
            encodeSeqToBits(currentSeqStr, encodedBitSequence);

            chunk[id] = std::make_tuple(encodedBitSequence, taxa->seqLen());
            if (_substitutionChunks->size() == 0){
                _rowIndexToId->push_back(id);
                _rowIndexToName->push_back(seqName);
            }

        }

        _substitutionChunks->push_back(chunk);

    }


    std::string getSequence(size_t id) {
        std::string fullSeq = "";
        for (auto &chunk: *_substitutionChunks) {
            std::array<char, MAX_COMPRESS_SIZE> encodedSeq = std::get<0>(chunk.at(id));
            size_t compressedSize = std::get<1>(chunk.at(id));
            size_t seqSize = std::get<2>(chunk.at(id));

            char decodedSeq[CHUNK_SIZE+1] = {0};
            decodedSeq[seqSize] = 0;

            size_t decompressedSize = HUF_PUBLIC_API::HUF_decompress(decodedSeq, seqSize,
                                                        encodedSeq.data(), compressedSize);



            fullSeq += decodedSeq;

        }

        return fullSeq;
    }

    size_t getIdFromRowIndex(size_t rowIndex) {
        return (*_rowIndexToId)[rowIndex];
    }

    std::string getNameFromRowIndex(size_t rowIndex) {
        return (*_rowIndexToName)[rowIndex];
    }

    bool isEmpty(){
        return (_substitutionChunks->size() == 0);
    }

    void encodeSeqToBits(const std::string &seq, std::vector<ALPHACHAR> &charVector) {
        if (_alphabet->size() == 4) {

        }
    }


    ~substitutionContainer(){}
};

#endif