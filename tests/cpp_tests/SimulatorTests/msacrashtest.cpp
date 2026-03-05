#include "../../../src/IndelSimulator.h"
#include "../../../src/SubstitutionSimulator.h"
#include "../../../libs/pcg/pcg_random.hpp"

// takes 10 minutes currently
int main() {
    tree tree_("../../trees/normalbranches_nLeaves10000.treefile");

    std::time_t t1 = 100;//std::time(0);

    SimulationContext<pcg64_fast> simContext(&tree_, t1);
    
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);

    // Zipf distribution with 'a' paramater of 1.7
    DiscreteDistribution d1({0.5095431711405822, 0.1568303071269451, 0.07871814521803239, 0.048270189115626465, 0.03303178006434935, 0.024228351178513602, 0.01864293638146846, 0.014856893415201527, 0.012160984068723036, 0.010166722872265421, 0.008645963958677104, 0.007457149799496436, 0.006508445711724147, 0.0057380367435973965, 0.005103003253080922, 0.004572745331946558, 0.004124941742034806, 0.003742981899246224, 0.003414285738097972, 0.0031291760165538885, 0.0028801040941869123, 0.0026611075564267036, 0.0024674233163113463, 0.00229520707878155, 0.0021413268904716906, 0.002003209144386965, 0.0018787222832818682, 0.0017660879701119473, 0.0016638125244592776, 0.0015706334865778777, 0.0014854775957006798, 0.0014074274672687627, 0.0013356949616762574, 0.0012695997452698497, 0.0012085519131340143, 0.0011520378136270332, 0.0010996084148699372, 0.0010508697030056355, 0.0010054747147966972, 0.0009631168927883048, 0.0009235245168354022, 0.0008864560163525095, 0.0008516960068952035, 0.0008190519253511694, 0.0007883511621300632, 0.000759438607808425, 0.0007321745468481026, 0.0007064328431279026, 0.0006820993717676363, 0.0006590706595873658});

    
    // d1.setSeed(t1);
    fill(insertionDists.begin(), insertionDists.end(), &d1);
    fill(deletionDists.begin(), deletionDists.end(), &d1);

    vector<double> insertionRates(tree_.getNodesNum() - 1);
    vector<double> deletionRates(tree_.getNodesNum() - 1);

    fill(insertionRates.begin(), insertionRates.end(), 0.03);
    fill(deletionRates.begin(), deletionRates.end(), 0.09);

    SimulationProtocol protocol(simContext.getTree()->getNodesNum() - 1);
    simContext.setProtocol(&protocol);

    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);
    protocol.setMaxInsertionLength(50);
    protocol.setMinSequenceSize(0);
    protocol.setSiteRateModel(SiteRateModel::INDEL_AWARE);


    int rootLength = 30000;
    protocol.setSequenceSize(rootLength);

    IndelSimulator<pcg64_fast> indelSim(simContext, &protocol);


    auto eventMap = indelSim.generateSimulation();
    std::cout << "done with event generation\n";

    modelFactory mFac;
    mFac.setReplacementModel(modelCode::NUCJC);
    // mFac.setModelParameters({0.25,0.25,0.25,0.25,0.1,0.2,0.3,0.4,0.5,0.6});

    mFac.setSiteRateModel({0.3068391864731156, 1.6931608135268843},
                          {0.5, 0.5},
                          {
                            {0.5, 0.5},
                            {0.5, 0.5}
                          });


    if (!mFac.isModelValid()) return 1;
    mFac.buildReplacementModel(); // to force model building
    
    auto categorySampler = mFac.getRateCategorySampler(protocol.getMaxInsertionLength());
    simContext.setCategorySampler(&categorySampler);

    auto msa = MSA<pcg64_fast>(eventMap, simContext);
    std::cout <<  msa.getMSAlength() << "\n";
    
    // std::cout << msa << "\n";
    



    return 0;


}