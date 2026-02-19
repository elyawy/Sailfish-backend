#include <chrono>

#include "../../../src/IndelSimulator.h"
#include "../../../src/SubstitutionSimulator.h"
#include "../../../libs/pcg/pcg_random.hpp"


// takes 10 minutes currently
int main() {
    tree tree_("../../trees/normalbranches_nLeaves10.treefile");

    std::time_t t1 = 42;//std::time(0);

    // time general simulation setup

    auto start = std::chrono::high_resolution_clock::now();
    SimulationContext<pcg64_fast> simContext(&tree_, t1);
    
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);

    // Zipf distribution with 'a' paramater of 1.7
    DiscreteDistribution d1({0.500597013622622,0.15407680416466002,0.07733607404815325,0.04742269916833997,0.03245183409648796,0.023802968878603163,0.018315618393686197,0.014596047786709592,0.011947471092366838,0.009988223562671806,0.008494164935846923,0.007326223038978798,0.00639417555007358,0.005637292815586038,0.005013408742738382,0.004492460672380736,0.004052519265065782,0.0036772655722418293,0.0033543403993036617,0.0030742364096058624,0.00282953749580224,0.002614385926699555,0.0024241022418638207,0.0022549096413393907,0.0021037311601300934,0.001968038376604936,0.0018457371577210268,0.001735080389074546,0.001634600615111441,0.00154305754135923,0.001459396750517745,0.0013827169647433223,0.0013122438819642994,0.0012473091132895959,0.0011873331108893934,0.0011318112414127814,0.0010803023567287186,0.0010324193607649647,0.0009878213819912308,0.000946207245246948,0.0009073100010352715,0.000870892320076723,0.0008367425994772478,0.0008046716569959321,0.0007745099135856162,0.0007461049831118826,0.0007193196030519133,0.000694029851883145,0.0006701236084401137,0.000647499216243446,0.0006260643220715613,0.0006057348631513914,0.0005864341815220788,0.0005680922475578417,0.0005506449774670669,0.0005340336319283628,0.0005182042849716617,0.0005031073538361537,0.0004886971818952107,0.00047493166787853947,0.0004617719355816668,0.00044918203906341153,0.00043712869901847034,0.00042558106659524523,0.0004145105114255328,0.0004038904310565982,0.0003936960793390335,0.00038390441163517177,0.00037449394498066375,0.0003654446315627292,0.0003567377440781168,0.00034835577170658,0.00034028232558560196,0.0003325020528024772,0.00032500055803344507,0.000317764332058729,0.0003107806864690648,0.0003040376939552896,0.00029752413363926596,0.00029122944096306576,0.0002851436617049981,0.00027925740973663626,0.00027356182817527154,0.0002680485536218626,0.00026270968320613465,0.0002575377441885291,0.00025252566589363025,0.00024766675377188466,0.00024295466540621684,0.00023838338829779629,0.00023394721928100186,0.00022964074543174937,0.0002254588263460113,0.00022139657767671146,0.00021744935582738057,0.00021361274371013342,0.00020988253748378816,0.00020625473419539554,0.00020272552025516334,0.0001992912606808326});

    
    // d1.setSeed(t1);
    fill(insertionDists.begin(), insertionDists.end(), &d1);
    fill(deletionDists.begin(), deletionDists.end(), &d1);

    vector<double> insertionRates(tree_.getNodesNum() - 1);
    vector<double> deletionRates(tree_.getNodesNum() - 1);

    fill(insertionRates.begin(), insertionRates.end(), 1.0);
    fill(deletionRates.begin(), deletionRates.end(), 1.0);

    SimulationProtocol protocol(simContext.getTree()->getNodesNum() - 1);
    simContext.setProtocol(&protocol);

    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);
    protocol.setMaxInsertionLength(150);
    protocol.setMinSequenceSize(10);
    protocol.setIndelRateModel(SiteRateModel::INDEL_AWARE);


    int rootLength = 10;
    protocol.setSequenceSize(rootLength);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "General simulation setup took " << duration << " microseconds.\n";

    std::cout << "Starting indel simulation...\n";
    // time indel simulator initialization in microseconds
    start = std::chrono::high_resolution_clock::now();
    IndelSimulator<pcg64_fast> indelSim(simContext, &protocol);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Indel simulator initialization took " << duration << " microseconds.\n";

    //time event generation in microseconds

    start = std::chrono::high_resolution_clock::now();
    auto eventMap = indelSim.generateSimulation();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Indel simulation took " << duration << " microseconds.\n";

    modelFactory mFac;
    // time model setup in microseconds
    start = std::chrono::high_resolution_clock::now();

    mFac.setReplacementModel(modelCode::LG);
    // mFac.setModelParameters({0.25,0.25,0.25,0.25,0.1,0.2,0.3,0.4,0.5,0.6});
    mFac.setSiteRateModel({0.01, 0.0,2.0,4.0},
                          {0.25, 0.25, 0.25, 0.25},
                          {
                            {0.59636944,0.27093557,0.11095153,0.02174347},
                            {0.2725427 ,0.35467176,0.26895975,0.10382579},
                            {0.11021473,0.27076059,0.35111379,0.2679109},
                            {0.0213646 ,0.10750572,0.2688647, 0.60226497}
                          });

    if (!mFac.isModelValid()) return 1;
    mFac.buildReplacementModel(); // to force model building
    end =  std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Model setup took " << duration << " microseconds.\n";

    // time substitution simulator initialization in microseconds
    start = std::chrono::high_resolution_clock::now();
    SubstitutionSimulator<pcg64_fast, 20> substitutionSim(mFac, simContext);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Substitution simulator initialization took " << duration << " microseconds.\n";
    
    //time MSA construction in microseconds
    start = std::chrono::high_resolution_clock::now();
    auto msa = MSA<pcg64_fast>(eventMap, simContext);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "MSA construction took " << duration << " microseconds.\n";
    std::cout << "MSA built. Number of sequences: " << msa.getNumberOfSequences() 
              << ", MSA length: " << msa.getMSAlength() << "\n";

    substitutionSim.setAlignedSequenceMap(msa);

    if (protocol.getSiteRateModel() == SiteRateModel::INDEL_AWARE) {
      auto rateCategories = msa.getPerSiteRateCategories();
      substitutionSim.setPerSiteRateCategories(rateCategories);
    }
    // substitutionSim.simulateAndWriteSubstitutions(msa.getMSAlength(), "output_newflowtest.fasta");
    // time substitution simulation in microseconds
    start = std::chrono::high_resolution_clock::now();
    auto fullContainer = substitutionSim.simulateSubstitutions(msa.getMSAlength());
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Substitution simulation took " << duration << " microseconds.\n";
    msa.fillSubstitutions(fullContainer);

    auto rateCategories = substitutionSim.getPerSiteRateCategories();

    // Print rate categories for each sites above the characters in the MSA
    std::cout << "Rate for all sites: ";
    for (auto& category: *rateCategories) {
        std::cout << category << " ";
    }
    std::cout << "\n";
    



    // start = std::chrono::high_resolution_clock::now();
    // msa.printFullMsa();
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    // std::cout << "MSA printing took " << duration << " microseconds.\n";


    return 0;


}