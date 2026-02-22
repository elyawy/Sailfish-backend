#include <iostream>
#include <cassert>


#include "../../../src/AvlTreeWithRates.h"

const size_t MAX_INSERTION_LENGTH = 100;

int main() {
    // simple test to check the way events are mapped to the tree structure.
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(9, 0);
    parentRates[0] = 1;
    // parentRates[1] = 1;
    parentRates[8] = 1;
    // print the parent rates
    std::cout << "Parent rates: ";
    for (size_t i = 0; i < parentRates.size(); ++i) {
        std::cout << parentRates[i] << " ";
    }
    std::cout << "\n";
    tree.init_tree(10, parentRates);
    // create a sampler with 4 actual categories:
    // such that if you are at 0 you transition to 1 with 100%.
    // if you are at 1 you transition to 2 with 100%, etc.
    CategorySampler sampler({{0.001, 0.999},
                             {0.999, 0.001}},
                            {0.5, 0.5}, MAX_INSERTION_LENGTH);

    std::mt19937_64 rng(42);

     std::vector<Event> events = {
        {DELETION, 7, 2},
        {INSERTION, 6, 1},
    };

    for (auto& ev : events) {
        std::cout << "\nEvent: " << (ev.type == INSERTION ? "INSERTION" : "DELETION")
                << " at position " << ev.position << ", size " << ev.length << "\n";
        bool success = tree.handle_event(ev, sampler, rng);
        std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
        std::cout << "Tree after event:\n" << tree.print_avl() << "\n";
    }

    
    std::vector<size_t> evolvedRateCategories;
    evolvedRateCategories.resize(tree.getTotalLength()-1);

    size_t processedSites = 0;
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        
        size_t length = it.val().length;
        size_t insertion = it.val().insertion;
        size_t siteInParent = it.key();
        if (insertion > 0 && siteInParent == 0 && length == 1 ) {
            // push vector at end of new vector
            std::cout << "adding insertions:\n";
            for (size_t i = 0; i < insertion; ++i) {
                evolvedRateCategories[processedSites] = it.val().rateCategories[i]; // new insertions get category 0
                processedSites++;
                std::cout << processedSites << " ";

            }
            std::cout << "\n";
            continue;
        }
        if (siteInParent == 0){
            siteInParent = 1;
            length -= 1;
        } // special case for the first site, which is at position 0 but corresponds to the first site in the parent sequence.
        
        std::cout << "copying from parent:\n";
        for (size_t pos_in_block = 0; pos_in_block < length; ++pos_in_block) {
            evolvedRateCategories[processedSites] = parentRates[(siteInParent-1) + pos_in_block];
            processedSites++;
            std::cout << "position in parent: " << siteInParent + pos_in_block << ", processed sites: " << processedSites << "\n";
        }
        std::cout << "\n";

        std::cout << "adding insertions:\n";
        for (size_t i = 0; i < insertion; ++i) {
            evolvedRateCategories[processedSites] = it.val().rateCategories[i]; // new insertions get category 0
            processedSites++;
            std::cout << processedSites << " ";
        }
        std::cout << "\n";
    }
    // there should be 14 rates categories in the evolved sequence,
    // corresponding to the 9 actual original sites and the 5 insertions.
    std::cout << "Evolved rate categories: ";
    for (size_t cat : evolvedRateCategories) {
        std::cout << cat << " ";
    }
    std::cout << "\n";


    return 0;
}