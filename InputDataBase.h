//
// Created by doby on 21/09/17.
//

#ifndef ADJUSTPVALUE_INPUTDATABASE_H
#define ADJUSTPVALUE_INPUTDATABASE_H

#include "InputData.h"
#include <string>
#include <vector>
#include "mysql_connection.h"
#include <cppconn/driver.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>
#include "PValue.h"

using namespace std;

class InputDataBase: public InputData {
    sql::Driver *driver;
    sql::Connection *con;
    sql::Statement *stmt;
    sql::ResultSet *res;
    sql::PreparedStatement *pstmt;
    vector<string> top_hits_ids;
public:
    InputDataBase (string user, string password, string database_name);
    ~InputDataBase();
    vector<vector<unsigned short>> createGenotypeMatrix(int lowInd, int upInd);
    vector<string> createPhenotypeStringVector();
    vector<unsigned short> createPhenotypeShortVector();
    vector<TableEntry> createTopHitsVector();
    void writeAdjustedPValuesToDB(vector<double> pvalues);
    void writeTopHitsTableToDB(vector<TableEntry> const & TopHitsTable);
};


#endif //ADJUSTPVALUE_INPUTDATABASE_H
