//
// Created by doby on 21/09/17.
//

#include "InputDataBase.h"
#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iterator>

InputDataBase::InputDataBase (string user, string password, string database_name){
    // Create a connection
    this->driver = get_driver_instance();
    this->con = driver->connect("tcp://127.0.0.1:3306", user, password);
    // Connect to the MySQL test database
    this->con->setSchema(database_name);
    this->stmt = this->con->createStatement();
}

InputDataBase::~InputDataBase(){
    delete this->res;
    delete this->stmt;
    delete this->con;
}

// Read data from table genotype: genotype strings for all SNPs in given interval and
// turn genotype strings into a vector<vector<unsigned short>>
// genotype will be used as input data for PValue::adjustPValue(tests, genotype, phenotype, parameters);
vector<vector<unsigned short>> InputDataBase::createGenotypeMatrix(int lowInd, int upInd){
    string genotype_string;
    vector <unsigned short> genotype_vector;
    vector<vector<unsigned short>> genotype;

    //lowInd--; upInd--; // Numeration starts with '0' in DataBase, numeration starts with '1' for user
    string statement = "SELECT snp_genotype FROM genotype2 limit " + to_string(lowInd) + ", " + to_string(upInd-lowInd+1);
    this->res = this->stmt->executeQuery(statement);
    while (this->res->next()) {
        genotype_string = this->res->getString("snp_genotype");
        genotype_vector = {};
        genotype_vector.reserve(genotype_string.size()); //to save on memory reallocations
        transform(begin(genotype_string), end(genotype_string), back_inserter(genotype_vector),
                       [](char c) {return c - '0';}); // Transform string to vector of unsigned short
        genotype.push_back(genotype_vector);
    }

    /*
    // Write resulting matrix to console
    for (vector<vector<unsigned char>>::size_type i = 0; i < genotype.size(); i++) {
        for (vector<unsigned char>::size_type j = 0; j < genotype[i].size(); j++ ) {
            cout << genotype[i][j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
    */

    return genotype;
}

// Read data from table phenotype to vector of strings
// phenotype will be used as input data for PValue::adjustPValue(tests, genotype, phenotype, parameters);
vector<string> InputDataBase::createPhenotypeStringVector(){
    this->res = this->stmt->executeQuery("SELECT phenotype_string FROM phenotype");
    vector<string> phen_vector;
    this->res->next();
    stringstream ss;
    ss.str(this->res->getString("phenotype_string"));
    string feature;
    while (getline(ss, feature, ',')) {
        phen_vector.push_back(feature);
    }

    /*
    // Write to console
    for (vector<string>::iterator it = phen_vector.begin(); it != phen_vector.end(); ++it) {
        cout << *it << endl;
    }
    cout << '\n';
     */

    return phen_vector;
}

vector<unsigned short> InputDataBase::createPhenotypeShortVector(){
    this->res = this->stmt->executeQuery("SELECT phenotype_string FROM phenotype2");
    vector<unsigned short> phen_vector;
    this->res->next();
    string phenotype_string = this->res->getString("phenotype_string");
    phen_vector.reserve(phenotype_string.size()); //to save on memory reallocations
    transform(begin(phenotype_string), end(phenotype_string), back_inserter(phen_vector),
              [](char c) {return c - '0';});
    /*
    // Write to console
    for (vector<unsigned char>::iterator it = phen_vector.begin(); it != phen_vector.end(); ++it) {
        cout << *it << endl;
    }
    cout << '\n';
     */

    return phen_vector;
}

// Read data from table density: test ids, indices of SNPs in the genotype matrix
// into a vector of TestsData structures
// tests will be used as input data for PValue::adjustPValue(tests, genotype, phenotype, parameters);
vector<TestsData> InputDataBase::createTopHitsVector(){
    vector<TestsData> tests;
    this->top_hits_ids = {};
    this->res = this->stmt->executeQuery("SELECT snp_id, test_id, lower, upper FROM density");
    while (this->res->next()) {
        tests.push_back({this->res->getString("test_id"), this->res->getInt("lower"), this->res->getInt("upper")});
        this->top_hits_ids.push_back(this->res->getString("snp_id"));
    }

    /*
    // Write to console
    for (vector<TestsData>::iterator it = tests.begin(); it != tests.end(); ++it) {
        cout << ' ' << (*it).ID << (*it).lower << (*it).upper << endl;
    }
    cout << '\n';
    */

    return tests;
}

void InputDataBase::writeAdjustedPValuesToDB(vector<double> pvalues){
    this->pstmt = this->con->prepareStatement("UPDATE density SET adjusted_pv = ? WHERE snp_id = ?");
    for (vector<double>::size_type j = 0; j < pvalues.size(); j++){
        this->pstmt->setDouble(1, pvalues[j]);
        this->pstmt->setString(2, this->top_hits_ids[j]);
        this->pstmt->executeUpdate();
    }
}

// Creates new table TopHits in DB and writes all the results
void InputDataBase::writeTopHitsTableToDB(vector<TableEntry> TopHitsTable){
    this->res = this->stmt->executeQuery("CREATE TABLE IF NOT EXISTS TopHits (test_index VARCHAR(20), "
                                                      "chromosome_index VARCHAR(20), ID VARCHAR(20), lower INT, "
                                                      "upper INT, adjusted_p_value DOUBLE)");
    this->pstmt = this->con->prepareStatement("INSERT INTO TopHits (test_index,chromosome_index,ID,lower,upper,"
                                                      "adjusted_p_value) VALUES(?,?,?,?,?,?)");
    for (vector<TableEntry>::size_type j = 0; j < TopHitsTable.size(); j++){
        this->pstmt->setString(1, TopHitsTable[j].test_index);
        this->pstmt->setString(2, TopHitsTable[j].chromosome_index);
        this->pstmt->setString(3, TopHitsTable[j].ID);
        this->pstmt->setInt(4, TopHitsTable[j].lower);
        this->pstmt->setInt(5, TopHitsTable[j].upper);
        this->pstmt->setDouble(6, TopHitsTable[j].adjusted_p_value);
        this->pstmt->executeUpdate();
    }
}