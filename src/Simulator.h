
#include <stack>
#include <random>

#include "../libs/Phylolib/includes/stochasticProcess.h"

#include "SimulationProtocol.h"
#include "BlockTree.h"
#include "MSA.h"
#include "Sequence.h"
#include "rateMatrixSim.h"
#include "modelFactory.h"


class Simulator
{
private:
    SimulationProtocol* _protocol;
    std::unique_ptr<rateMatrixSim> _substitutionSim;
    size_t _seed;
	std::mt19937_64 _mt_rand;
    std::uniform_real_distribution<double> _biased_coin;
    std::shared_ptr<std::vector<bool>> _nodesToSave;
    // std::uniform_int_distribution<int> _fair_die;

public:
    Simulator(SimulationProtocol* protocol): _protocol(protocol),
    _seed(protocol->getSeed()), _mt_rand(protocol->getSeed()),
    _biased_coin(0,1) {
        // std::cout << "simulator ready!\n";
        DiscreteDistribution::setSeed(_seed);
        _nodesToSave = std::make_shared<std::vector<bool>>(_protocol->getTree()->getNodesNum(), false);
        setSaveStateLeaves(_protocol->getTree()->getRoot());
    }

    void initSimulator() {
        _seed = _protocol->getSeed();
        DiscreteDistribution::setSeed(_seed);
        _mt_rand.seed(_seed);
    }

    void resetSimulator(SimulationProtocol* newProtocol) {
        _protocol = newProtocol;
        initSimulator();
    }

    std::vector<BlockMap> runSimulator(size_t numberOfSimulation) {
        // std::vector<MSA*> msas;
        std::vector<BlockMap> msaPrecursors;

        for (size_t i = 0; i < numberOfSimulation; ++i) {
            BlockMap blockmap = generateSimulation();
            msaPrecursors.push_back(blockmap);
        }
        return msaPrecursors;
    }

    BlockMap generateSimulation() {
        size_t sequenceSize = _protocol->getSequenceSize();
        tree::TreeNode *rootNode = _protocol->getTree()->getRoot();

        BlockMap nodeToBlockMap;
        BlockList rootBlockList;
        std::array<size_t, 3> rootBlock = {0,sequenceSize + 1, 0};
        rootBlockList.push_back(rootBlock);
        nodeToBlockMap[rootNode->id()] = std::make_tuple(rootBlockList, sequenceSize + 1);
        generateIndelsRecursively(nodeToBlockMap, *rootNode);
        return nodeToBlockMap;
    }

    void generateIndelsRecursively(BlockMap &blockmap,const tree::TreeNode &currentNode) {
        if (currentNode.isLeaf()) {
            return;
        }
        auto blockTuple = blockmap.at(currentNode.id());
        size_t seqLength = std::get<static_cast<int>(BLOCKLIST::LENGTH)>(blockTuple);
        size_t correctedSeqLength = seqLength - 1;

        for (size_t i = 0; i < currentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode =  currentNode.getSon(i);
            auto newBlockTuple = simulateAlongBranch(correctedSeqLength, childNode->dis2father(), childNode->id()-1);
            blockmap[childNode->id()] = newBlockTuple;
            generateIndelsRecursively(blockmap, *(currentNode.getSon(i)));
        }
    }


    std::tuple<BlockList, size_t> simulateAlongBranch(size_t seqSize, double branchLength, size_t nodePosition) {
        size_t sequenceSize = seqSize;
        size_t minSequenceSize = _protocol->getMinSequenceSize();
        // std::cout << "sequenceSize=" << sequenceSize << "\n";

        BlockTree blocks(sequenceSize);

        double insertionRate = _protocol->getInsertionRate(nodePosition);
        double deletionRate = _protocol->getDeletionRate(nodePosition);
        // std::cout << "insertionRate=" << insertionRate << "\n";
        // std::cout << "deletionRate=" << deletionRate << "\n";

        DiscreteDistribution* insertionLengthDistribution = _protocol->getInsertionDistribution(nodePosition);
        DiscreteDistribution* deletionLengthDistribution = _protocol->getDeletionDistribution(nodePosition);

        double sampledDeletionLength = deletionLengthDistribution->drawSample();

        double sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
        double sequenceWiseDeletionRate = 1.0 * deletionRate * (sequenceSize + (sampledDeletionLength - 1));
        if (sequenceSize <= minSequenceSize) sequenceWiseDeletionRate = 0.0;
        // std::cout << "sequenceWiseInsertionRate=" << sequenceWiseInsertionRate << "\n";
        // std::cout << "sequenceWiseDeletionRate=" << sequenceWiseDeletionRate << "\n";

        double lambdaParam = sequenceWiseInsertionRate + sequenceWiseDeletionRate;
        // std::cout << "lambda=" << lambdaParam << "\n";
        std::exponential_distribution<double> distribution(lambdaParam);


        double waitingTime = distribution(_mt_rand);
        // std::cout << "waitingTime=" << waitingTime << "\n";
        // std::cout << "branchLength=" << branchLength << "\n";

        while (waitingTime < branchLength) {

            double insertionProbability = sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate);

            int eventIndex;
            size_t eventLength;
            event eventType;

            double coinFlip = _biased_coin(_mt_rand);

            if (coinFlip < insertionProbability) {
                auto _fair_die = std::uniform_int_distribution<int>(0, sequenceSize);
                eventIndex = _fair_die(_mt_rand);
                // std::cout << eventIndex << " ";
                eventLength = insertionLengthDistribution->drawSample();
                eventType = event::INSERTION;
            } else {
                auto _fair_die = std::uniform_int_distribution<int>(1 - (sampledDeletionLength-1), sequenceSize);
                eventIndex = _fair_die(_mt_rand);
                eventLength = sampledDeletionLength;
                if (eventIndex < 1) {
                    eventLength = eventLength + (eventIndex-1);
                    eventIndex = 1;
                }
                // std::cout << eventIndex << " ";

                if (eventLength + eventIndex > sequenceSize) {
                    eventLength = sequenceSize - eventIndex + 1;
                }
                eventType = event::DELETION;
            }
            blocks.handleEvent(eventType, eventIndex, eventLength);

            // if (!blocks.checkLength()){
            //     std::cout << blocks.printTree() << "\n";
            //     std::cout << "eventType=" << eventType
            //               << " eventIndex=" << eventIndex
            //               << " sequenceSize=" << sequenceSize
            //               << " eventLength=" << eventLength << "\n";
            //     errorMsg::reportError("length error!");
            // }


            sequenceSize = blocks.length() - 1;
            sampledDeletionLength = deletionLengthDistribution->drawSample();

            branchLength = branchLength - waitingTime;
            sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
            sequenceWiseDeletionRate =  1.0 * deletionRate * (sequenceSize + (sampledDeletionLength - 1));
            if (sequenceSize <= minSequenceSize) sequenceWiseDeletionRate = 0.0;

            lambdaParam = sequenceWiseInsertionRate + sequenceWiseDeletionRate;
            std::exponential_distribution<double> distribution(lambdaParam);
            waitingTime = distribution(_mt_rand);

        }
        return std::make_tuple(blocks.getBlockList(), blocks.length());
    }

    void initSubstitionSim(modelFactory& mFac) {
        _substitutionSim = std::make_unique<rateMatrixSim>(mFac, _nodesToSave);
        // _substitutionSim->setSeed(_seed);
        _substitutionSim->setRng(&_mt_rand);
    }

    std::vector<double> getSiteRates() {
        std::vector<double> temp = _substitutionSim->getSiteRates();
        _substitutionSim->clearRatesVec();
        return temp;
    }


    void setNodesToSave(std::vector<size_t> nodeIDs) {
        std::fill(_nodesToSave->begin(), _nodesToSave->end(), false);
        for(auto &nodeID: nodeIDs) {
            (*_nodesToSave)[nodeID] = true;
        }
    }

    void setSaveAllNodes() {
        for (size_t i = 0; i < _nodesToSave->size(); i++) {
            (*_nodesToSave)[i] = true;
        }
    }

    void setSaveRoot() {
        (*_nodesToSave)[0] = true;
    }


    void changeNodeSaveState(size_t nodeID) {
        (*_nodesToSave)[nodeID] = !(*_nodesToSave)[nodeID];
    }

    bool getNodeSaveState(size_t nodeID) {
        return (*_nodesToSave)[nodeID];
    }

    const std::vector<bool> getNodesSaveList() {
        return (*_nodesToSave);
    }

    void setSaveStateLeaves(const tree::nodeP &node) {
        for(auto &node: node->getSons()) {
            if (node->isLeaf()) (*_nodesToSave)[node->id()] = true;
            setSaveStateLeaves(node);
        }
    }

    void setSaveRates(bool saveRates) {
        _substitutionSim->setSaveRates(saveRates);
    }

    void writeSeqContainerToFiles(const std::unique_ptr<sequenceContainer> seqContainer, const std::string& directory, bool isFirst) {
        
        auto seq = seqContainer->constTaxaBegin();

        for(; seq != seqContainer->constTaxaEnd(); ++seq) {

            std::ofstream file(directory + "/" + std::to_string((*seq).id()) + ".fasta", 
                               isFirst ? std::ios_base::out : std::ios_base::app);
            if (!file.is_open()) {
                // Handle error - couldn't open the file
                continue;
            }
            if (isFirst) file << ">" << (*seq).name() << "\n";
            file << (*seq).toString();
            file.close();

        }
    }

    void simulateAndWriteSubstitutions(size_t sequenceLength, const std::string& directory) {
        size_t chunkSize = 1024;
        size_t numberOfChunks = (size_t)sequenceLength / chunkSize;
        size_t remainder = sequenceLength % chunkSize; /* Likely uses the result of the division. */
        // std::cout << "number of chunks: " << numberOfChunks << ", with remainder: " << remainder << "\n"; 

        for (size_t i = 0; i < numberOfChunks; i++) {
            // std::cout << i << "\n";
            _substitutionSim->generate_substitution_log(chunkSize);
            auto seqContainer = _substitutionSim->getSequenceContainer();
            if (i == 0) {
                writeSeqContainerToFiles(std::move(seqContainer), directory, true);
                continue;
            }
            writeSeqContainerToFiles(std::move(seqContainer), directory, false);
        }
        if (remainder == 0) return;
        
        _substitutionSim->generate_substitution_log(remainder);
        auto seqContainer = _substitutionSim->getSequenceContainer();

        if (numberOfChunks == 0) {
            writeSeqContainerToFiles(std::move(seqContainer), directory, true);
            return;
        }
        writeSeqContainerToFiles(std::move(seqContainer), directory, false);
    }

    std::shared_ptr<sequenceContainer> simulateSubstitutions(size_t sequenceLength) {
        size_t chunkSize = 1024;
        size_t numberOfChunks = (size_t)sequenceLength / chunkSize;
        size_t remainder = sequenceLength % chunkSize; /* Likely uses the result of the division. */
        // std::cout << "number of chunks: " << numberOfChunks << ", with remainder: " << remainder << "\n"; 

        std::shared_ptr<sequenceContainer> fullSequence;
        for (size_t i = 0; i < numberOfChunks; i++) {
            // std::cout << i << "\n";
            _substitutionSim->generate_substitution_log(chunkSize);
            auto seqContainer = _substitutionSim->getSequenceContainer();
            if (i == 0) {
                fullSequence = std::move(seqContainer);
                continue;
            }
            for (size_t j = 0; j < fullSequence->numberOfSeqs(); ++j) {
                int idOfSeq = seqContainer->placeToId(j);
                (*fullSequence)[idOfSeq] += (*seqContainer)[idOfSeq];
            }
        }
        if (remainder == 0) return fullSequence;

        
        _substitutionSim->generate_substitution_log(remainder);
        auto seqContainer = _substitutionSim->getSequenceContainer();

        if (fullSequence == nullptr) {
            fullSequence = std::move(seqContainer);
            return fullSequence;
        }
        for (size_t j = 0; j < fullSequence->numberOfSeqs(); ++j) {
            int idOfSeq = seqContainer->placeToId(j);
            (*fullSequence)[idOfSeq] += (*seqContainer)[idOfSeq];
        }
        

        return fullSequence;

    }


    ~Simulator(){}
};
