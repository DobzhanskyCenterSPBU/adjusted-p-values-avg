#include <iostream>
#include "PValue.h"
#include "InputFile.h"
#include <chrono>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "PerformanceTests.h"
#include <boost/math/special_functions/gamma.hpp>
#include "InputDataBase.h"
#include <plinkio/plinkio.h>
#include "InputPLINK.h"
#include "AuxiliaryFunctions.cpp"
//#include "omp.h"

using namespace std;
using namespace std::chrono;

vector<TableEntry> readTopHitsTable(string filename, int number_of_lines) {
    ifstream handle;
    handle.open(filename);
    if (!handle.is_open()) {
        cerr << "No such file!" << endl;
        exit(1); // What exit code should be used? How to skip to the next iteration?
    }

    TableEntry currentEntry;
    vector<TableEntry> result;
    string data, rank, test_index, chr_index, snp_ind_left, snp_ind_center, snp_ind_right, lower_pos, center_pos,
            upper_pos, count, neg_log_dens, dens_p_value, density;

    getline(handle, data); // Skip the first line
    for (int i = 0; i < number_of_lines; ++i) // What other way is there? while(handle.good()){}
    {
        getline(handle, rank, ';');
        getline(handle, test_index, ';'); //!
        getline(handle, chr_index, ';'); //!
        getline(handle, snp_ind_left, ';');
        getline(handle, snp_ind_center, ';');
        getline(handle, snp_ind_right, ';');
        getline(handle, lower_pos, ';'); //!
        getline(handle, center_pos, ';');
        getline(handle, upper_pos, ';'); //!
        getline(handle, count, ';');
        getline(handle, neg_log_dens, ';');
        getline(handle, dens_p_value, ';'); //!
        getline(handle, density);

        // double naive_p_value = boost::lexical_cast<double>(dens_p_value); // ???????????????????

        currentEntry = {.test_index = test_index, .chromosome_index = chr_index, .ID = "", .lower = stoi(lower_pos),
                .upper = stoi(upper_pos), .adjusted_p_value = -2};
        result.push_back(currentEntry); // Is int enough for index ?????
    }

    return result;
}

void insertTestIDInTopHitsTable (vector<TableEntry> & top_hits, vector<string> IDs) {


    if (top_hits.size() != IDs.size()) {
        cerr << "The array of IDs is not the same length as the Top Hits Table!" << endl;
        exit(2);
    }

    for (typename vector<string>::size_type i = 0; i < IDs.size(); i++) {
        top_hits[i].ID = IDs[i];
    }

}

void insertAdjPValInTopHitsTable (vector<TableEntry> & top_hits, vector<double> adj_val) {


    if (top_hits.size() != adj_val.size()) {
        cerr << "The array of adjusted p values is not the same length as the Top Hits Table!" << endl;
        exit(2);
    }

    for (typename vector<string>::size_type i = 0; i < adj_val.size(); i++) {
        top_hits[i].adjusted_p_value = adj_val[i];
    }

}

int main(int argc, char* argv[]) {

    // Run unit tests
    testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    //simplePerformanceTest(10000000, "cd");
    //runMultiplePerformanceTests(10);
    //runMultiplePerformanceTestsMultipleTimes(100);

    /*
     * Manual data initialization + genotype is read from .txt
     * Results are presented in standard output
     */
    /*
    vector<TableEntry>  tests = {{"","","cd",2,5,0}, {"","","r",6,8,0}, {"","","d",10,11,0}, {"","","a",12,14,0}, {"","","d",16,19,0}};
    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0};
    InputFile genotype("test.txt");
    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = false;

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Print results
    cout << endl << "Adjusted P-values:" << endl;
    printVector(results);
     */



    /*
     * Input from DataBase
     * Results are placed into a database
     */
    /*
    // Data initialization
    InputDataBase DBConnectObj(argv[1], argv[2], argv[3]);

    vector<TableEntry> tests = DBConnectObj.createTopHitsVector(); // Genotype

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = false;

    vector<unsigned short> phenotype = DBConnectObj.createPhenotypeShortVector();
    //vector<string> phenotypeVector = DBConnectObj.createPhenotypeStringVector(); // Use if phenotype is in string format in file
    //vector<unsigned short> phenotype = PValue::mapPhenotypeValuesToChar(phenotypeVector);

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, DBConnectObj, phenotype, parameters);

    // Print results
    cout << endl << "Adjusted P-values:" << endl;
    printVector(results);

    // Write results to database
    DBConnectObj.writeAdjustedPValuesToDB(results);
    */



    /*
     * Input data from .csv and PLINK files
     * Results are placed into a database
     */

    char* PLINK_filename = "/Users/doby/Desktop/plink_file_work/input/cleaned_all_final";
    string topHitsTable_filename = "F_density_50000_test.csv";
    int number_of_lines = 10;
    vector<string> ID = {"cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd"};

    vector<TableEntry> top_hits_table;
    InputPLINK PLINK_data(PLINK_filename); // Genotype object
    vector <unsigned short> phenotype;
    ExecutionParameters parameters; //Will be passed through the command line
    InputDataBase DBConnectObj("root", "root1234", "test_top_hits_output");

    top_hits_table = readTopHitsTable(topHitsTable_filename, number_of_lines);
    insertTestIDInTopHitsTable(top_hits_table, ID);
    phenotype = PLINK_data.createPhenotypeVector();
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = false;

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(top_hits_table, PLINK_data, phenotype, parameters);

    // Write results into DB
    insertAdjPValInTopHitsTable(top_hits_table, results); // Insert adjusted P Values into Top Hits Table
    DBConnectObj.writeTopHitsTableToDB(top_hits_table);


    return 0;
}