
#include <stack>
#include <random>
#include <memory>

#include "SimulationContext.h"
#include "SimulationProtocol.h"
#include "Event.h"

template<typename RngType = std::mt19937_64>
class IndelSimulator
{
private:
    SimulationProtocol* _protocol;
    tree* _tree;
    tree::TreeNode* _rootNode;
    std::uniform_real_distribution<double> _biased_coin;
    const std::vector<bool>& _nodesToSave;
    RngType& _rng;

public:
    IndelSimulator(SimulationContext<RngType>& simContext, SimulationProtocol* protocol):
    _protocol(protocol),
    _tree(simContext.getTree()),
    _rootNode(simContext.getRoot()),
    _biased_coin(0,1),
    _nodesToSave(simContext.getNodesToSave()),
    _rng(simContext.getRng()) {}


    void updateSimulationProtocol(SimulationProtocol* newProtocol) {
        _protocol = newProtocol;
    }


    EventMap generateSimulation() {
        size_t sequenceSize = _protocol->getSequenceSize();
        tree::TreeNode *rootNode = _rootNode;

        EventMap nodeToEventMap;
        nodeToEventMap.resize(_tree->getNodesNum());
        // set root sequence size by inserting a dummy event list with a single insertion of length sequenceSize
        nodeToEventMap[rootNode->id()] = EventSequence{{event::INSERTION, 0, sequenceSize}};
        generateIndelsRecursively(nodeToEventMap, sequenceSize, *rootNode);
        return nodeToEventMap;
    }

    void generateIndelsRecursively(EventMap &eventmap, size_t parentSequenceLength ,const tree::TreeNode &currentNode) {
        if (currentNode.isLeaf()) {
            return;
        }

        for (size_t i = 0; i < currentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode =  currentNode.getSon(i);
            auto eventTuple = simulateAlongBranch(parentSequenceLength, childNode->dis2father(), childNode->id()-1);
            // get the new sequence length after indels
            size_t newSequenceLength = eventTuple.second;
            // store the events for the current node
            eventmap[childNode->id()] = eventTuple.first;
            generateIndelsRecursively(eventmap, newSequenceLength, *(currentNode.getSon(i)));
        }
    }


    std::pair<EventSequence, size_t> simulateAlongBranch(size_t seqSize, double branchLength, size_t nodePosition) {
        EventSequence events;
        
        size_t sequenceSize = seqSize;
        size_t minSequenceSize = _protocol->getMinSequenceSize();

        double insertionRate = _protocol->getInsertionRate(nodePosition);
        double deletionRate = _protocol->getDeletionRate(nodePosition);

        DiscreteDistribution* insertionLengthDistribution = _protocol->getInsertionDistribution(nodePosition);
        DiscreteDistribution* deletionLengthDistribution = _protocol->getDeletionDistribution(nodePosition);

        double sampledDeletionLength = deletionLengthDistribution->drawSample(_rng);

        double sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
        double sequenceWiseDeletionRate = 1.0 * deletionRate * (sequenceSize + (sampledDeletionLength - 1));
        if (sequenceSize <= minSequenceSize) sequenceWiseDeletionRate = 0.0;

        double lambdaParam = sequenceWiseInsertionRate + sequenceWiseDeletionRate;
        std::exponential_distribution<double> distribution(lambdaParam);


        double waitingTime = distribution(_rng);

        while (waitingTime < branchLength) {

            double insertionProbability = sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate);

            int eventIndex;
            size_t eventLength;
            event eventType;

            double coinFlip = _biased_coin(_rng);

            if (coinFlip < insertionProbability) {
                auto _fair_die = std::uniform_int_distribution<int>(0, sequenceSize);
                eventIndex = _fair_die(_rng);

                eventLength = insertionLengthDistribution->drawSample(_rng);
                eventType = event::INSERTION;
            } else {
                auto _fair_die = std::uniform_int_distribution<int>(1 - (sampledDeletionLength-1), sequenceSize);
                eventIndex = _fair_die(_rng);
                eventLength = sampledDeletionLength;
                if (eventIndex < 1) {
                    eventLength = eventLength + (eventIndex-1);
                    eventIndex = 1;
                }

                if (eventLength + eventIndex > sequenceSize) {
                    eventLength = sequenceSize - eventIndex + 1;
                }
                eventType = event::DELETION;
            }
            events.push_back({eventType, static_cast<size_t>(eventIndex), eventLength});
            // std::cout << "  Event: " << (eventType == INSERTION ? "Insertion" : "Deletion") << " at " << eventIndex << " length " << eventLength << "\n";

            if (eventType == INSERTION) {
                sequenceSize += eventLength;
            } else {  // DELETION
                sequenceSize -= eventLength;
            }

            sampledDeletionLength = deletionLengthDistribution->drawSample(_rng);

            branchLength = branchLength - waitingTime;
            sequenceWiseInsertionRate = 1.0 * insertionRate * (sequenceSize + 1);
            sequenceWiseDeletionRate =  1.0 * deletionRate * (sequenceSize + (sampledDeletionLength - 1));
            if (sequenceSize <= minSequenceSize) sequenceWiseDeletionRate = 0.0;

            lambdaParam = sequenceWiseInsertionRate + sequenceWiseDeletionRate;
            std::exponential_distribution<double> distribution(lambdaParam);
            waitingTime = distribution(_rng);

        }

        return std::make_pair(events, sequenceSize);
    }

    ~IndelSimulator(){}
};
