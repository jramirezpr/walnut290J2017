#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <random>     // For test case
#include <functional> // For test case
using namespace std;

#include "SimDb_SQL.h"

// This sample program demonstrates the use of SimDb_SQL class
// to create and use an SQLite database for the specification
// of input and the capturing of output for the execution
// of a program or routine.
//
// In this example, the routine being run is one that evaluates
// an approximation to the integral of a function F(x) over [-1,1] using
// a random sampling. The parameters in the routine being varied
// are the number of samples used and the seed for the
// random number generator.
//
// Step 1: (a) Connect to an SQLite database (creating the SQLite database
//         file if necessary).
//
//         (b) Defining the database structure that will hold simulation
//         input and output data.
//
//         (c) Store input data values in the database by inserting
//         a row in a database table for each instance the program or
//         routine is to be run.
//
//         (d) Disconnect from the database. This step is not
//         necessary if the program that creates the database
//         is also used to execute the program or routine with
//         the specified data inputs.
//
// Step 2: (a) Connect to an existing database that contains the input
//             parameters (if needed).
//
//         (b) For each row in the database extract input parameters and
//             run the simmulation or routine. Capture the output and then set
//             the output data values in the database.
//
//         (c) Close the database
//
//
// Step 3: Using another program work with the database that is produced.
//
//         Other than having a general understanding of what SQL database tables,
//         columns, and rows are, you don't need to know much SQL to use the SimDb_SQL
//         class. To work with the databases produced, you will need to know
//         some SQL. A good way to learn is to install the firefox browser
//         add-on SQLite manager from
//
//         https://addons.mozilla.org/en-US/firefox/addon/sqlite-manager/
//
//         use it to open up the database created with this sample, and experiment
//         with browsing the database and executing SQL statements.
//
//
// Created for Math 290J : 11/21/2015 Chris Anderson

// Forward declaration of the routine being tested.
// The routine source is at the end of main.


double randomIntegralApprox(double a, double b, std::function<double(double)> F, long randomSeed, long sampleCount);


int main(int argc, char* argv[])
{
    SimDb_SQL simDb;
    simDb.setVerboseFlag();

// simDb.setDBverboseFlag(); // Will display each SQL command being executed

//################################################################
//   Step 1 (a),(b):  Create database and specify structure.
//################################################################

    // Connects to an sqlite file "First.db".
    // If the file doesn't exist it will be created.

    simDb.connectDataBase("First.db");

    // Specify a table to hold the input and output data values.

    string tableName = "RunData";

    // Set up the structure of the table by specifying the column names
    // and data type affinity (e.g. the data for the values in the column).

    vector <string> columnNames;
    vector <string> columnDataTypes;

    columnNames.push_back("status");            columnDataTypes.push_back("TEXT");
    columnNames.push_back("Sample_Size");       columnDataTypes.push_back("INTEGER");
    columnNames.push_back("Random_Seed");       columnDataTypes.push_back("INTEGER");
    columnNames.push_back("Integral_Est");      columnDataTypes.push_back("REAL");

    // Create the table. If the table already exists then this member function
    // will append to the existing table the columns specfied by columnNames that
    // are not already present in the table. Existing table columns are not
    // altered.
    //

    simDb.createTable(tableName,columnNames,columnDataTypes);

    // Clear any existing rows (table data) in the table because we want to
    // specify a new set of input data. If the table is a new table,
    // then this is a no-op.

    simDb.clearTableRows(tableName);

//################################################################
//   Step 1 (c):  Specify input data
//################################################################

    vector <long > seeds       = {12345,314159,0};
    vector <long > sampleCount = {100,200,400,800,1600,3200};

    long SimId;

    for(int i = 0; (int) i < seeds.size(); i++)
    {
    for(int j = 0; (int) j < sampleCount.size(); j++ )
    {
    	SimId = simDb.insertRow(tableName); // Capture the SimId so it can be used for updating values
    	simDb.updateRow(SimId,tableName,"status","none");
    	simDb.updateRow(SimId,tableName,"Sample_Size",sampleCount[j]);
    	simDb.updateRow(SimId,tableName,"Random_Seed",seeds[i]);
    }}

//################################################################
//   Step 1 (d):  Close the database.
//################################################################

    // This is just being done so that in Step 2, a non-trivial
    // connection statement can be demonstrated.

    simDb.disconnectDataBase();

//################################################################
//   Step 2 (a): Connect to the database
//################################################################

	// Connecting to an existing database. Specifying noCreate so
	// that if this database doesn't exist, and error message
	// will be generated.

	bool noCreate = true;
    simDb.connectDataBase("First.db",noCreate);

//################################################################
//   Step 2 (b): Loop over rows in the database, extracting
//               parameters, running the routine, capturing
//               output and then updating the database.
//################################################################

    vector <long > SimId_Values = simDb.getSimIdValues(tableName);
    long               seedData;
    long         sampleSizeData;
    double     integralEstimate;

    // Specify the function to be tested. In this case it's (3/2)x^2
    // which has the integral of 1 over [-1,1].

    double a = -1.0;
    double b =  1.0;

	std::function<double(double)> F = [](double x) {return (3.0/2.0)*x*x;};

    for(long i = 0; i < SimId_Values.size(); i++)
    {
    SimId = SimId_Values[i];
    printf("Processing Task : %-10ld \n",SimId); 

    seedData         = simDb.getIntegerData(SimId,tableName,"Random_Seed");
    sampleSizeData   = simDb.getIntegerData(SimId,tableName,"Sample_Size");

    // Run the routine can capture the output

    integralEstimate = randomIntegralApprox(a,b,F, seedData,sampleSizeData);

    // Update the database

    simDb.updateRow(SimId,tableName,"status","done");
    simDb.updateRow(SimId,tableName,"Integral_Est",integralEstimate);
    }


//################################################################
//   Step 2 (c): Disconnect the database
//################################################################

    simDb.disconnectDataBase();

    printf("XXX Execution Completed XXXX\n");

	return 0;
}



//################################################################
//
//   Routine that uses a monte-carlo estimate of the
//   value of the integral of F(x) over [a,b]
//
//################################################################

double randomIntegralApprox(double a, double b, std::function<double(double)> F, long randomSeed, long sampleCount)
{
	mt19937_64                           randomGenerator;
	uniform_real_distribution<double>      distribution;
	randomGenerator.seed(randomSeed);

	// Initialize the distribution to be uniform in the interval [-1,1]

	uniform_real_distribution<double>::param_type distParams(-1.0,1.0);
	distribution.param(distParams);

    double xVal;
    double dx = (b-a)/(double)sampleCount;
    double integralApprox = 0.0;

	for(long i = 0; i < sampleCount; i++)
	{
	    xVal = distribution(randomGenerator);
	    integralApprox += F(xVal)*dx;
    }

    return integralApprox;
}

