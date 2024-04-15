#include <ctime>

#include "../../src/Simulator.h"
// #include "definitions.h"



// takes 10 minutes currently
int main() {
    tree tree_("../trees/normalbranches_nLeaves10000.treefile");
    // tree tree_("(A:0.1,B:0.2);", false);

    // tree_.getRoot()->orderSonsByHeight();
    std::time_t t1 = 12;//std::time(0);
    // DiscreteDistribution::setSeed(t1);
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);

    // Zipf distribution with 'a' paramater of 1.7
    DiscreteDistribution d1({0.500597013622622,0.15407680416466002,0.07733607404815325,0.04742269916833997,0.03245183409648796,0.023802968878603163,0.018315618393686197,0.014596047786709592,0.011947471092366838,0.009988223562671806,0.008494164935846923,0.007326223038978798,0.00639417555007358,0.005637292815586038,0.005013408742738382,0.004492460672380736,0.004052519265065782,0.0036772655722418293,0.0033543403993036617,0.0030742364096058624,0.00282953749580224,0.002614385926699555,0.0024241022418638207,0.0022549096413393907,0.0021037311601300934,0.001968038376604936,0.0018457371577210268,0.001735080389074546,0.001634600615111441,0.00154305754135923,0.001459396750517745,0.0013827169647433223,0.0013122438819642994,0.0012473091132895959,0.0011873331108893934,0.0011318112414127814,0.0010803023567287186,0.0010324193607649647,0.0009878213819912308,0.000946207245246948,0.0009073100010352715,0.000870892320076723,0.0008367425994772478,0.0008046716569959321,0.0007745099135856162,0.0007461049831118826,0.0007193196030519133,0.000694029851883145,0.0006701236084401137,0.000647499216243446,0.0006260643220715613,0.0006057348631513914,0.0005864341815220788,0.0005680922475578417,0.0005506449774670669,0.0005340336319283628,0.0005182042849716617,0.0005031073538361537,0.0004886971818952107,0.00047493166787853947,0.0004617719355816668,0.00044918203906341153,0.00043712869901847034,0.00042558106659524523,0.0004145105114255328,0.0004038904310565982,0.0003936960793390335,0.00038390441163517177,0.00037449394498066375,0.0003654446315627292,0.0003567377440781168,0.00034835577170658,0.00034028232558560196,0.0003325020528024772,0.00032500055803344507,0.000317764332058729,0.0003107806864690648,0.0003040376939552896,0.00029752413363926596,0.00029122944096306576,0.0002851436617049981,0.00027925740973663626,0.00027356182817527154,0.0002680485536218626,0.00026270968320613465,0.0002575377441885291,0.00025252566589363025,0.00024766675377188466,0.00024295466540621684,0.00023838338829779629,0.00023394721928100186,0.00022964074543174937,0.0002254588263460113,0.00022139657767671146,0.00021744935582738057,0.00021361274371013342,0.00020988253748378816,0.00020625473419539554,0.00020272552025516334,0.0001992912606808326});

    
    // d1.setSeed(t1);
    fill(insertionDists.begin(), insertionDists.end(), &d1);
    fill(deletionDists.begin(), deletionDists.end(), &d1);

    vector<double> insertionRates(tree_.getNodesNum() - 1);
    vector<double> deletionRates(tree_.getNodesNum() - 1);

    fill(insertionRates.begin(), insertionRates.end(), 0.03);
    fill(deletionRates.begin(), deletionRates.end(), 0.09);
    // fill(insertionRates.begin(), insertionRates.end(), 0.01);
    // fill(deletionRates.begin(), deletionRates.end(), 0.01);

    SimulationProtocol protocol(&tree_);


    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);

    int rootLength = 8000;
    protocol.setSequenceSize(rootLength);

    protocol.setSaveAncestral(false);

    protocol.setSeed(t1);

    Simulator sim(&protocol);

    // sim.initSimulator();
    using BlockMap = std::map<std::string, BlockTree>;


    std::vector<BlockMap> blockmaps = sim.runSimulator(1);

    std::cout << "finished all indel simulations\n";
    std::vector<MSA> msas = MSA::generateMSAs(blockmaps, tree_.getRoot());
    int msaLength = msas[0].getMSAlength();

    std::cout << "length of the MSA will be: " << msaLength << "\n";
    // std::cin.get();

    modelFactory mFac(&tree_);

    mFac.setAlphabet(alphabetCode::NUCLEOTIDE);
    mFac.setReplacementModel(modelCode::NUCJC);
    // mFac.setModelParameters({0.25,0.25,0.25,0.25,0.1});

    mFac.setGammaParameters(1.0, 1); // TODO: ALLOW 1 CATEGORY!

    if (!mFac.isModelValid()) return 0;

    sim.initSubstitionSim(mFac);
    std::cout << "initializing subs sim" << "\n";
    // std::cin.get();
    

    // std::cout << "number of nodes to simulate: " << tree_.getNodesNum() - 1 << "\n";
    // auto seqContainer = sim.simulateSubstitutions(msaLength);

    auto fullContainer = sim.simulateSubstitutions(msaLength);
    // for (size_t i = 0; i < tree_.getLeavesNum(); i++)
    // {
    //    std::cout << i << "\n";
    //    size_t nodeId = (*fullContainer).placeToId(i);
    //    std::cout << (*fullContainer)[nodeId] << "\n";
    // }
    
    
    std::cout << "finished all substitutions" << "\n";
    // std::cin.get();


    msas[0].fillSubstitutions(fullContainer);
    std::cout << "filled MSA" << "\n";
    // std::cin.get();

    // msas[0].printFullMsa();
    // msas[0].writeFullMsa("/home/elyalab/fasta.fasta");

    // blockmaps.clear();
    // msas.clear();
    // std::cout << "finished generating MSAs\n";

    // std::cin.get();
    // std::vector<tree::nodeP> leaves;

    // tree_.getAllLeaves(leaves, tree_.getRoot());
    // for (auto &node: leaves)
    // {
    //     std::cout << blockmaps[0][node->name()].printTree();
    // }

    return 0;


}