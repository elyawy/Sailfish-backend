
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
    // std::uniform_int_distribution<int> _fair_die;

    public:
    using BlockMap = std::map<std::string, BlockTree>;
    Simulator(SimulationProtocol* protocol): _protocol(protocol),
    _seed(protocol->getSeed()), _mt_rand(protocol->getSeed()),
    _biased_coin(0,1) {
        // std::cout << "simulator ready!\n";
        DiscreteDistribution::setSeed(_seed);
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

    std::map<std::string, BlockTree> generateSimulation() {
        using nodeP = tree::nodeP;

        size_t sequenceSize = _protocol->getSequenceSize();
        std::stack<nodeP> nodes;
        nodes.push(_protocol->getTree()->getRoot());
        nodeP currentNode = nodes.top();

        BlockMap nodeToBlockMap;
        nodeToBlockMap[currentNode->name()] = BlockTree(sequenceSize);
        size_t nodePosition = 0;
        while (!nodes.empty()) {
            nodes.pop();
            if (!currentNode->isLeaf()) {
                for (auto node: currentNode->getSons()) {
                    nodes.push(node);
                }
            } else {
                if (nodes.empty()) break;
            }
            currentNode = nodes.top();
            sequenceSize = nodeToBlockMap[currentNode->father()->name()].length() - 1;
            BlockTree blocks = simulateAlongBranch(sequenceSize, currentNode->dis2father(), nodePosition);
            nodeToBlockMap[currentNode->name()] = blocks;

            ++nodePosition;
        }

        return nodeToBlockMap;
    }

    BlockTree simulateAlongBranch(size_t seqSize, double branchLength, size_t nodePosition) {
        size_t sequenceSize = seqSize;

        // std::cout << "sequenceSize=" << sequenceSize << "\n";

        BlockTree blocks(sequenceSize);

        double insertionRate = _protocol->getInsertionRate(nodePosition);
        double deletionRate = _protocol->getDeletionRate(nodePosition);

        // std::cout << "insertionRate=" << insertionRate << "\n";
        // std::cout << "deletionRate=" << deletionRate << "\n";

        DiscreteDistribution* insertionLengthDistribution = _protocol->getInsertionDistribution(nodePosition);
        DiscreteDistribution* deletionLengthDistribution = _protocol->getDeletionDistribution(nodePosition);

        double sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
        double sequenceWiseDeletionRate = 1.0 * deletionRate * sequenceSize;

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

            size_t eventIndex;
            size_t eventLength;
            event eventType;

            double coinFlip = _biased_coin(_mt_rand);

            if (coinFlip < insertionProbability) {
                auto _fair_die = std::uniform_int_distribution<int>(1, sequenceSize + 1);
                eventIndex = _fair_die(_mt_rand);
                // std::cout << eventIndex << " ";
                eventLength = insertionLengthDistribution->drawSample();
                eventType = event::INSERTION;
            } else {
                auto _fair_die = std::uniform_int_distribution<int>(1, sequenceSize);
                eventIndex = _fair_die(_mt_rand);
                // std::cout << eventIndex << " ";

                eventLength = deletionLengthDistribution->drawSample();
                if (eventLength + eventIndex > sequenceSize) {
                    eventLength = sequenceSize - eventIndex;
                }
                eventType = event::DELETION;
            }
            // std::cout << "event type " << eventType << "\tposition " << eventIndex << "\tlength " << eventLength << "\n";
            blocks.handleEvent(eventType, eventIndex, eventLength);                


            sequenceSize = blocks.length() - 1;

            branchLength = branchLength - waitingTime;
            sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
            sequenceWiseDeletionRate = 1.0 * deletionRate * sequenceSize;
            
            lambdaParam = sequenceWiseInsertionRate + sequenceWiseDeletionRate;
            std::exponential_distribution<double> distribution(lambdaParam);
            waitingTime = distribution(_mt_rand);

        }
        
        return blocks;
    }

    void initSubstitionSim(modelFactory& mFac) {
        _substitutionSim = std::make_unique<rateMatrixSim>(mFac);
        // _substitutionSim->setSeed(_seed);
        _substitutionSim->setRng(&_mt_rand);
    }

    std::shared_ptr<sequenceContainer> simulateSubstitutions(size_t sequenceLength) {
        _substitutionSim->generate_substitution_log(sequenceLength);


        std::shared_ptr<sequenceContainer> sharedSeqContainer = _substitutionSim->toSeqDataWithoutInternalNodes();
        _substitutionSim->initSim();
        return sharedSeqContainer;
    }


    ~Simulator(){}
};
