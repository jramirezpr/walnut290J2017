#include <iostream>
#include <string>
#include <cstdio>
#include <map>
#include <vector>
using namespace std;

#include "sqlite3.h"
#include "SimDb_Utilities.h"

/*
Class SimDb_SQL

Beta Version : 15.11.20

This is a C++ interface class to a SQLite database file and provides a
core set of member functions that are useful for opening, closing, storing data,
and retrieving data in SQLite database files.

The SimDb_SQL class was created to eliminate some of the coding complexity
associated with using the C interface to SQLite databases that is part of the
standard SQLite distribution obtainable from https://www.sqlite.org.

The class member functions allow one to perform various database operations
with arguments that are standard C++ data type and data structures; thus
circumventing the need to wrestle with the problems associated with
data type conversion necessary when working the the SQLite C interface.

Specification of the conversions internally used in this class:

(*) Double and float data are stored with 15 digits of precision
    as REAL data in the SQLite database file.

(*) Long and int data are both stored as INTEGER data
    in an SQLite database file.

(*) stl::string data is stored as TEXT in the SQLite database file.

Error checking (mostly type checking) is performed using a cached version of the
database schema so that checks can be performed without issuing
underlying database calls.

The database connection statement in this class opens the database
in serialzed mode (multi-thread support, no restrictions). This requires
that the sqlite3 library that is linked to the executable be compiled with
multi-threaded support. See https://www.sqlite.org/threadsafe.html

After one has created a simulation databases, one typically works with
the data in the database using other programs. Generally, the use of these
programs require some familiarity with SQL.

SQLite command line program    : sqlite3

SQLite Manager firefox plug-in : https://addons.mozilla.org/en-US/firefox/addon/sqlite-manager/

Python Pandas                  : http://pandas.pydata.org/

Veusz with the sqlite3 plugin  : https://github.com/waveform80/veusz_plugins

================================================================================

Creation Date: Nov. 20, 2015
Author       : Chris Anderson

 */

 /*
#############################################################################
#
# Copyright  2015 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

#ifndef _SimDb_SQL_
#define _SimDb_SQL_

class SimDb_SQL
{
public :

    SimDb_SQL()
    {
    connectionFlag          = false;
	db                      = 0;
	existingDBflag          = false;
	verboseFlag             = false;
	verboseDBflag           = false;
	haltOnErrorFlag         = true;
	defaultInitFlag         = true;

    }

	virtual ~SimDb_SQL()
	{
	disconnectDataBase();
	}

    void setVerboseFlag()     {verboseFlag     = true;}
    void clearVerboseFlag()   {verboseFlag     = false;}
    void setDBverboseFlag()   {verboseDBflag   = true;}
    void clearDBverboseFlag() {verboseDBflag   = false;}
    void setHaltOnErrors()    {haltOnErrorFlag = true;}
    void clearHaltOnErrors()  {haltOnErrorFlag = false;}

//  If the SQLite database file exists, this routine connects to the existing
//  file. Rows inserted are appended to the database file.

//  If the SQLite database file doesn't exist, and noCreate = false (the default),
//  then an SQLite database file with the specified name is created.

	int connectDataBase(const string& databaseFileName, bool noCreate = false)
	{
	string dataBase;
	existingDBflag = sqlUtil.fileExists(databaseFileName.c_str());

	if((noCreate)&&(not existingDBflag))
	{
    if(verboseFlag)
	{
	printf("Database file does not exist : %s\n",databaseFileName.c_str());
	}
	return 1;
	}
	connectionFlag = false;

	// Open and/or open and create the database.
	// Using sqlite3_open_v2 to be explicit about opening the
	// database in serialized mode
	//
	// https://www.sqlite.org/c3ref/open.html

    const char* nullFilePtr = 0;

    utf8chars.clear();
    sqlUtil.createUTF8characterEncoding(databaseFileName,utf8chars);
	int rc = sqlite3_open_v2(&utf8chars[0], &db,SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE|SQLITE_OPEN_FULLMUTEX,nullFilePtr);

    if(rc)
	{
	printf("SimDb_SQL Error: \n");
	printf("SQ_ERROR: %s \n",sqlite3_errmsg(db));
	if(db) {sqlite3_close(db);}
	return 1;
	}

	if(verboseFlag)
	{
	printf("Successfully Opened SQLite3 Database File : %s\n",databaseFileName.c_str());
	}

    dataBaseName.assign(databaseFileName);
	connectionFlag = true;

	// Cache the current schema

	updateDbSchemaData();

    return 0;
	}

    // This member function disconnects from an attached database

	void disconnectDataBase()
	{
	if(not connectionFlag) return;


    if(db) {sqlite3_close(db);}
    if(verboseFlag)
	{
	printf("Successfully Closed SQLite3 Database \n");
	}
    connectionFlag = false;
    existingDBflag = false;
    }
    //
    // Creates a database table with a SimId column and columns specified by columnNames
    // with data types columnDataTypes.
    //
    // Admissible data type specification : INTEGER, REAL or TEXT
    //
	int createTable(string& tableName,vector < string > columnNames, vector < string > columnDataTypes)
	{
    nameCheck(tableName);

	if(not connectionFlag) return 1;
    ostringstream sqCmdStream;
    vector < char >     sqCmd;

    vector < string > existingColumnNames;
    vector < string > existingColumnDataTypes;
    long updateCount;
    string columnType;

    long  columnSize = columnNames.size();
    if(columnSize != (long)columnDataTypes.size())
    {
    printf("SimDb_SQL Error: \n");
	printf("Create Table : Insufficient column data type sepecification \n");
	printf("Column names specified : \n");
	for(int i = 0; i < columnSize; i++)                  {printf("%d : %s\n",i,columnNames[i].c_str());}
	printf("Column data type specified : \n");
	for(int i = 0; i < (int)columnDataTypes.size(); i++) {printf("%d : %s\n",i,columnDataTypes[i].c_str());}
	disconnectDataBase();
	return 1;
	}

    for(long k = 0; k < columnSize; k++)
    {
    	nameCheck(columnNames[k]);
        dataTypeCheck(columnDataTypes[k]);
    }

    sqCmdStream.str("");

    if(not isTable(tableName.c_str()))
    {
        sqCmdStream.str("");
    	sqCmdStream << "CREATE TABLE " << tableName.c_str() << "(";

    	sqCmdStream << "SimId" << " "  << "INTEGER PRIMARY KEY";
    	if(columnSize > 0) sqCmdStream << ", ";

    	for(long k = 0; k < columnSize; k++)
    	{
    		sqCmdStream << columnNames[k].c_str() << " " << columnDataTypes[k].c_str();
    		if(k < columnSize-1)
    		{
    		sqCmdStream << ", ";
    		}
    	}

    	sqCmdStream << ");";
        executeDBput(sqCmdStream.str());
    }
    else // Alter existing tables by adding any columns as needed (never drop a column)
    {
    	updateCount   = 0;
    	existingColumnNames.clear();
    	getDBcolumnInfo(tableName,existingColumnNames,existingColumnDataTypes);

    	for(long k = 0; k < columnSize; k++)
    	{
    		if(not sqlUtil.isInArray(existingColumnNames,columnNames[k]))
    		{
    		updateCount++;
    		}
    	}
    	if(updateCount != 0) //ALTER TABLE table_name ADD column_name datatype
    	{
    	for(long k = 0; k < columnSize; k++)
    	{
    		if(not sqlUtil.isInArray(existingColumnNames,columnNames[k]))
    		{
    		    sqCmdStream.str("");
    		    sqCmdStream << "ALTER TABLE " << tableName.c_str() << " ADD ";
    			sqCmdStream << columnNames[k].c_str() << " " << columnDataTypes[k].c_str();
    		    sqCmdStream << ";";
    		    executeDBput(sqCmdStream.str());
    		}
    	}
    	}
    }

   // Cache internal schema information

   updateDbSchemaData();

   return 0;
   }

   // Removes all rows from the specified table and resets
   // the SimId counter to 0, so the next insertion will
   // have SimId = 1.

   int clearTableRows(const string& tableName)
   {
   if(not isTable(tableName)) { tableError(tableName); return 1;}

   string sqCommandString;
   sqCommandString =  "DELETE FROM " + tableName + ";";
   executeDBput(sqCommandString);
   return 0;
   }

   // Removes the specified row.

   // If a row less than the maximal SimId is deleted, there is
   // no re-ordering of the SimId values. The next insertion will
   // have a SimId that is the current maximal SimId + 1.
   //
   // If the last row is deleted then the next insertion will
   // have the same SimId as the last row.
   //

   int deleteRow(long SimId, const string& tableName)
   {
    if(not isTable(tableName)) { tableError(tableName); return 1;}

    ostringstream sqCmdStream;
    vector < char >     sqCmd;

    sqCmdStream.str("");
    sqCmdStream << "DELETE FROM " << tableName.c_str() << " WHERE SimId = " << SimId << ";";
    executeDBput(sqCmdStream.str());
   }

    //
    // Inserts an empty row into the database and returns the SimId value
    // of the inserted row.
    //
    // Returns SimId = index of inserted row (automatically incremented)
    //
    // Returns -1 if an insertion error occurs.

    int insertRow(const string& tableName)
    {
    if(not isTable(tableName)) { tableError(tableName); return 1;}

    vector< string > colNames;
    vector< string >  colData;
    vector< string >  colType;
    int SimId            = -1;

    insertAutoIncrementedRow(tableName,colNames,colData,colType,SimId);
    return SimId;
    }

    int updateRow(long SimId, const string& tableName,
    const  string& name, const double& value)
    {
    updateRowDataCheck(tableName, name, "REAL");
    string valueString;
    valueString.assign(sqlUtil.getSQdataString(value).c_str());
    return updateRow(SimId, tableName, name, valueString);
    }

    int updateRow(long  SimId, const string& tableName,
    const  string& name, const float& value)
    {
    updateRowDataCheck(tableName, name, "REAL");
    string valueString;
    valueString.assign(sqlUtil.getSQdataString(value).c_str());
    return updateRow(SimId, tableName, name, valueString);
    }

    int updateRow(long SimId, const string& tableName,
    const  string& name, const int& value)
    {
    updateRowDataCheck(tableName, name, "INTEGER");
    string valueString;
    valueString.assign(sqlUtil.getSQdataString(value).c_str());
    return updateRow(SimId, tableName, name, valueString);
    }

    int updateRow(long SimId, const string& tableName,
    const  string& name, const long& value)
    {
    updateRowDataCheck(tableName, name, "INTEGER");
    string valueString;
    valueString.assign(sqlUtil.getSQdataString(value).c_str());
    return updateRow(SimId, tableName, name, valueString);
    }

    int updateRow(long SimId, const string tableName, const
    string colName, const string& colData)
	{
	if(not isTable(tableName))          { tableError(tableName); return 1;}
	if(not isColumn(tableName,colName)) { columnError(tableName,colName); return 1;}

	bool textFlag;
	ostringstream tmpStream;
	tmpStream << "UPDATE " << tableName.c_str() << " SET ";

    tmpStream << colName.c_str() << " = ";
    textFlag = false;
	if((getColumnDateType(tableName,colName)).compare("TEXT") == 0){textFlag = true;}

	if(textFlag) tmpStream << "\'";
    tmpStream << colData.c_str();
    if(textFlag) tmpStream << "\'";

	tmpStream << " WHERE SimId = " <<  SimId << ";";
	return executeDBput(tmpStream.str());
    }


    // SELECT colName FROM tableName WHERE SimId = XXX
    //
    // colName refers to a column that has TEXT affinity in the database.
    // The data is returned as a string in colData. The value returned by
    // a null column is an empty string.

    string getStringData(const long SimId, const string& tableName, const string& colName)
    {
        if(not isTable(tableName))          { tableError(tableName); return string();}
	    if(not isColumn(tableName,colName)) { columnError(tableName,colName); return string();}

        string colDataType = getColumnDateType(tableName, colName);
        if(colDataType.compare("TEXT") != 0)
        {
            printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	 	    printf(" SIMDb_SQL ERROR : getStringData() accessing %s column data \n",colDataType.c_str());
	 	    printf(" Table  name      : %s \n",tableName.c_str());
	 	    printf(" Column name      : %s \n",colName.c_str());
	 	    printf(" Column data type : %s \n\n",colDataType.c_str());
	 	    errorHandling();
	 		return string();
        }
        string colData;

        ostringstream tmpStream;

        int rc             = 0;
        sqlite3_stmt* stmt = 0;
    	tmpStream << "SELECT " << colName.c_str() <<  " FROM " << tableName.c_str() << " WHERE SimId = " << SimId;
        if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

        sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
        rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);

        if(rc != SQLITE_OK)
        {
	 	    printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	 	    sqlite3_finalize( stmt );
	 	    errorHandling();
	 		return string();
        }

        const char* sPtr;
        if( sqlite3_step( stmt ) == SQLITE_ROW )
        {
        sPtr = (const char*)sqlite3_column_text(stmt, 0 );
    	if(sPtr != 0){colData.assign(sPtr);}
        }

        sqlite3_finalize( stmt );

        return colData;
    }

    // SELECT colName FROM tableName WHERE SimId = XXX
    //
    // colName refers to a column that has INTEGER affinity in the database.
    // The data is returned as a string in colData. The value returned by
    // a null column is 0.
    //
    long getIntegerData(const long SimId, const string& tableName, const string& colName)
    {
        if(not isTable(tableName))          { tableError(tableName); return 0;}
	    if(not isColumn(tableName,colName)) { columnError(tableName,colName); return 0;}

        string colDataType = getColumnDateType(tableName, colName);
        if(colDataType.compare("INTEGER") != 0)
        {
            printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	 	    printf(" SIMDb_SQL ERROR : getIntegerData() accessing %s column data \n",colDataType.c_str());
	 	    printf(" Table  name      : %s \n",tableName.c_str());
	 	    printf(" Column name      : %s \n",colName.c_str());
	 	    printf(" Column data type : %s \n\n",colDataType.c_str());
	 	    errorHandling();
	 		return 0;
        }
        long colData = 0;
        ostringstream tmpStream;
        int rc             = 0;
        sqlite3_stmt* stmt = 0;
	    const char* sPtr;

    	tmpStream << "SELECT " << colName.c_str() <<  " FROM " << tableName.c_str() << " WHERE SimId = " << SimId;
    	if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

        sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);

        rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);

        if(rc != SQLITE_OK)
        {
	 		printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	 	    sqlite3_finalize( stmt );
	 	    errorHandling();
	 		return 0;
        }

        if( sqlite3_step( stmt ) == SQLITE_ROW )
        {
        sPtr = (const char*)sqlite3_column_text(stmt, 0 );
    	if(sPtr != 0){colData = atoi(sPtr);}
        }

        sqlite3_finalize( stmt );
        return colData;
    }

    // SELECT colName FROM tableName WHERE SimId = XXX
    //
    // colName refers to a column that has REAL affinity in the database.
    // The data is returned as a string in colData. The value returned by
    // a null column is 0.0
    //

    double getRealData(const long SimId, const string& tableName, const string& colName)
    {
        if(not isTable(tableName))          { tableError(tableName); return 0.0;}
	    if(not isColumn(tableName,colName)) { columnError(tableName,colName); return 0.0;}

        string colDataType = getColumnDateType(tableName, colName);
        if(colDataType.compare("REAL") != 0)
        {
            printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	 	    printf(" SIMDb_SQL ERROR : getRealData() accessing %s column data \n",colDataType.c_str());
	 	    printf(" Table  name      : %s \n",tableName.c_str());
	 	    printf(" Column name      : %s \n",colName.c_str());
	 	    printf(" Column data type : %s \n\n",colDataType.c_str());
	 	    errorHandling();
	 		return 0.0;
        }

        double colData = 0.0;
        ostringstream tmpStream;
        int rc             = 0;
        sqlite3_stmt* stmt = 0;
	    const char* sPtr;

       	tmpStream << "SELECT " << colName.c_str() <<  " FROM " << tableName.c_str() << " WHERE SimId = " << SimId;
        if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

        sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
        rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);

        if(rc != SQLITE_OK)
        {
	 		printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	 	    sqlite3_finalize( stmt );
	 	    errorHandling();
	 		return 0.0;
        }

        if( sqlite3_step( stmt ) == SQLITE_ROW )
        {
        sPtr = (const char*)sqlite3_column_text(stmt, 0 );
    	if(sPtr != 0){colData = atof(sPtr);}
        }

        sqlite3_finalize( stmt );
        return colData;
    }

    vector< long > getSimIdValues(const string& tableName)
    {
    if(not isTable(tableName))          { tableError(tableName); return vector<long>();}

    vector <string > colName(1,"SimId");
    vector <string > colSpec;
    vector <string > colSpecVal;
    vector <string > colSpecType;
    string orderBy;
    vector < vector< string> > colData;
    getDbData(tableName, colName, colSpec, colSpecVal, colSpecType,  orderBy, colData);

    vector <long > SimIdValues(colData[0].size(),0);
    for(int i = 0; i < (int) colData[0].size(); i++)
    {
    SimIdValues[i] = atoi(colData[0][i].c_str());
    }

    return SimIdValues;
    }

    int getMaxSimId(const string& tableName)
    {
    if(not isTable(tableName))          { tableError(tableName); return 1;}

    long SimId = -1;
    getDBintegerColumnMax(tableName,"SimId", SimId);
    return SimId;
    }

    int updateRow(long SimId, const string tableName, const vector <string>& colNames,
    const vector <string>& colData)
	{
	if(not isTable(tableName))          { tableError(tableName); return 1;}

	long varCount = colNames.size();
	bool textFlag;
	ostringstream tmpStream;

	tmpStream << "UPDATE " << tableName.c_str() << " SET ";

	for(long k = 0; k < varCount; k++)
	{
	if(not isColumn(tableName,colNames[k])) { columnError(tableName,colNames[k]); return 1;}

    tmpStream << colNames[k].c_str() << " = ";
    textFlag = false;
	if((getColumnDateType(tableName,colNames[k])).compare("TEXT") == 0){textFlag = true;}

	if(textFlag) tmpStream << "\'";
    tmpStream << colData[k].c_str();
    if(textFlag) tmpStream << "\'";
    if(k != varCount-1){tmpStream << ", ";}
	}
	tmpStream << " WHERE SimId = " <<  SimId << ";";
	return executeDBput(tmpStream.str());
    }


    // Uses the internally cached information about database schema

    bool isTable(const string& table)
    {
    nameCheck(table);

    std::map < string, map < string, string > >::iterator it;
    for(it =  dbSchemaData.begin(); it != dbSchemaData.end(); ++it)
    {
    	if(table.compare(it->first) == 0) return true;
    }
    return false;
    }

    // Uses the internally cached information about database schema

    bool isColumn(const string& tableName, const string& column)
    {
    nameCheck(tableName);
    nameCheck(column);

    map< string, string > colData = dbSchemaData[tableName];
    map< string, string >::iterator it;

    for(it =  colData.begin(); it != colData.end(); ++it)
    {
    	if(column.compare(it->first) == 0) return true;
    }
    return false;
    }

	// This routine gets the all the table names in the database

    int getDBtableNames(vector < string > & tableNames)
    {
    tableNames.clear();
    if(connectionFlag == false){return 1;}
    ostringstream tmpStream;

    tableNames.clear();
	sqlite3_stmt* stmt = 0;
    tmpStream << "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
    if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

    sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
    int rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);
	if(rc != SQLITE_OK)
	{
	printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
    sqlite3_finalize( stmt );
    errorHandling();
	return 1;
	}
    const char* sPtr;
    while(sqlite3_step( stmt ) == SQLITE_ROW)
    {
    sPtr = (const char*)sqlite3_column_text(stmt, 0 );
    if(sPtr != 0) {tableNames.push_back(sPtr);}
    }
    sqlite3_finalize( stmt );

    return 0;
    }

    // Uses the internally cached information about database schema

    string getColumnDateType(const string& tableName, const string& column)
    {
    if(not isTable(tableName))          { tableError(tableName);         return string();}
    if(not isColumn(tableName,column))  { columnError(tableName,column); return string();}

    map< string, string > colData = dbSchemaData[tableName];
    return colData[column];
    }

    //
    // Error checking routine to insure that the data used to update
    // the particular column value has the correct type.
    //

    int errorHandling()
	{
    if(haltOnErrorFlag)
	{
	printf("XXXXXXX     Program Terminated        XXXXXXX\n");
	exit(1);
    }
    return 1;
	}


    // This routine executes a standard SQL database command
    // specified by dbCommand. Any output from the execution of the
    // statement is discarded,

	int executeDBput(const string& dbCommand)
	{
	if(connectionFlag == false){return 1;}

	char*        sqLiteErrMsg;
	int          rc;
	if(verboseDBflag)
	{
	printf("%s\n",dbCommand.c_str());
	}

	sqlUtil.createUTF8characterEncoding(dbCommand,utf8chars);
	rc = sqlite3_exec(db,&utf8chars[0], 0, 0, &sqLiteErrMsg);
	if(rc != SQLITE_OK)
	{
    		printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
    		printf(" Offending command  : %s \n ",dbCommand.c_str());
    		sqlite3_free(sqLiteErrMsg);
            errorHandling();
    		return 1;
	}

	return 0;
	}


    int dataTypeCheck(const string& name)
    {
    nameCheck(name);
    if(name.compare("TEXT") == 0)    return 0;
    if(name.compare("REAL") == 0)    return 0;
    if(name.compare("INTEGER") == 0) return 0;
    printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	printf(" SIMDb_SQL ERROR           :  Unacceptable data type specified. \n");
	printf(" Acceptable types          :  TEXT, REAL, INTEGER \n");
    printf(" Offending type specified  : %s \n", name.c_str());
    errorHandling();
	return 1;
    }


    int tableError(const string& name)
    {
    printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	printf(" SIMDb_SQL ERROR : Attempt to access non-existent table. \n");
    printf(" Offending name  : %s \n", name.c_str());
    errorHandling();
	return 1;
    }

    int columnError(const string& tableName, const string& name)
    {
    printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	printf(" SIMDb_SQL ERROR : Attempt to access non-existent column. \n");
	printf(" Table           : %s \n", tableName.c_str());
    printf(" Offending name  : %s \n", name.c_str());
    errorHandling();
	return 1;
    }

    int nameCheck(const string& name)
    {
    string nameStar = sqlUtil.trim(name);
    if(nameStar.compare(name) != 0)
    {
    printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	printf(" SIMDb_SQL ERROR :  Unacceptable name. \n");
	printf(" Name has leading or trailing whitespace. \n");
    printf(" Offending name  : %s \n", name.c_str());
    errorHandling();
	return 1;
    }

    if(name.find(" ") != string::npos)
    {
	printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
    printf(" SIMDb_SQL ERROR :  Unacceptable name. \n");
	printf(" Internal spaces not allowed.  \n");
    printf(" Offending name  : %s \n", name.c_str());
	errorHandling();
	return 1;
    }
    }

	int updateRowDataCheck(const string& tableName, const string& colName, const string& updateDataType)
	{
	if(not isTable(tableName))          { tableError(tableName);               return 1;}
    if(not isColumn(tableName,colName)) { columnError(tableName,colName);      return 1;}

    string colDataType = getColumnDateType(tableName, colName);
    if(colDataType.compare(updateDataType) != 0)
        {
            printf(" XXXXXXXXXXX ERROR XXXXXXXXXXXX \n");
	 	    printf(" SIMDb_SQL ERROR : updateRow() parameters data type mis-match \n");
	 	    printf(" Table  name      : %s \n",tableName.c_str());
	 	    printf(" Column name      : %s \n",colName.c_str());
	 	    printf(" Column data type : %s \n\n",colDataType.c_str());
	 	    printf(" Input  data type : %s \n\n",updateDataType.c_str());
	 	    errorHandling();
	 		return 1;
        }
	return 0;
	}

    // This routine gets the all the column names in a particular table.
    // In SQLite3 there doesn't seem to be a primitive for it, so the
    // sql constructor for the table is parsed for the information.
    //
    // The reason for using the procedure here instead of using
    // the function sqlite3_column_name() is that the latter is only
    // applicable after a select statement --- which will not work if
    // the table happens to be empty.
    //

    int getDBcolumnInfo(const string& tableName,vector < string >& columnNames,vector < string >& columnDataTypes)
    {


    columnNames.clear();
    if(connectionFlag == false){return 1;}

    ostringstream tmpStream;

    string sqlStatement;
	sqlite3_stmt* stmt = 0;

    tmpStream << "SELECT sql FROM sqlite_master WHERE tbl_name = \'" << tableName.c_str() << "\';";
    if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

    sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
    int rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);
	if(rc != SQLITE_OK)
	{
	printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	sqlite3_finalize( stmt );
    errorHandling();
	return 1;
	}
	const char* sPtr;
    if(sqlite3_step( stmt ) == SQLITE_ROW)
    {
    sPtr = (const char*)sqlite3_column_text(stmt, 0 );
    if(sPtr != 0){sqlStatement.assign(sPtr);}
    }
    sqlite3_finalize( stmt );

    // Parse the sql constructor to extract the column names
    // this assumes that the SQLite3 result has the form
    // "CREATE TABLE TableName(var1 var1Type [Other Attributes], var2 var2Type [Other Attributes], var3 var3Type [Other Attributes], ... )"

    // strip () from argument to sql statement

    int s1 = sqlStatement.find_first_of("(",0);
    int s2 = sqlStatement.find_last_of(")",sqlStatement.size()-1);
    string sqlArg = sqlStatement.substr(s1+1, s2-s1-1);

    // Tokenize using commas to obtain the elements in the table construction

    vector < string > sqlParts;
    sqlUtil.createTokens(sqlArg,sqlParts,",");

    // Tokenize with spaces and select the first two tokens to
    // be the column name and the column data type

    vector < string > colInfo;
    for(int i = 0; i < (int)sqlParts.size(); i++)
    {
    colInfo.clear();
    sqlUtil.createTokens(sqlParts[i],colInfo," ");
    columnNames.push_back(colInfo[0]);
    columnDataTypes.push_back(colInfo[1]);
    }

    return 0;
	}


    int insertAutoIncrementedRow(const string& tableName, const vector< string >& colNames,
	const vector< string >& colData, const vector< string > &colType, int& SimId)
	{
	long varCount = colNames.size();
	bool textFlag;
	ostringstream tmpStream;

	tmpStream << "INSERT INTO " << tableName.c_str() << " ";

	if(varCount == 0)
	{
	tmpStream << "DEFAULT values;";
	}
	else
	{
    tmpStream << "(";
	for(long k = 0; k < varCount; k++)
	{
    tmpStream << colNames[k].c_str();
    if(k != varCount-1){tmpStream << ", ";}
	}
	tmpStream << " ) VALUES ( ";
    for(long k = 0; k < varCount; k++)
	{
	textFlag = false;
	if(colType[k].compare("TEXT") == 0){textFlag = true;}

	if(textFlag) tmpStream << "\'";
    tmpStream << colData[k].c_str();
    if(textFlag) tmpStream << "\'";

    if(k != varCount-1){tmpStream << ", ";}
	}
	tmpStream << " );";
	}

    // Execute the insertion and query to obtain the last row added

    string dbCommand = tmpStream.str();

	if(connectionFlag == false){return 1;}

	char*        sqLiteErrMsg;
	int          rc;

	if(verboseDBflag) {printf("%s\n",dbCommand.c_str());}

	sqlUtil.createUTF8characterEncoding(dbCommand,utf8chars);

	rc = sqlite3_exec(db,&utf8chars[0], 0, 0, &sqLiteErrMsg);
	if(rc != SQLITE_OK)
	{
    		printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
    		printf(" Offending command  : %s \n ",dbCommand.c_str());
    		sqlite3_free(sqLiteErrMsg);
            errorHandling();
    		return -1;
	}

	SimId = sqlite3_last_insert_rowid(db);

    return SimId;
	}


	// SELECT MAX(column_name) FROM table_name

    int getDBintegerColumnMax(const string& tableName, const string& columnName, long& maxVal)
    {
    if(not isTable(tableName))              { tableError(tableName);               return 1;}
    if(not isColumn(tableName,columnName))  { columnError(tableName,columnName);   return 1;}

    if(connectionFlag == false){return 1;}
    ostringstream tmpStream;

	sqlite3_stmt* stmt = 0;
    tmpStream << "SELECT MAX(" << columnName.c_str() << ") FROM " << tableName.c_str() << ";";
    if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}

    sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
    int rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);
	if(rc != SQLITE_OK)
	{
	printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	sqlite3_finalize( stmt );
    errorHandling();
	return 1;
	}

    if(sqlite3_step( stmt ) == SQLITE_ROW)
    {
    maxVal = sqlite3_column_int(stmt, 0 );
    }
    sqlite3_finalize( stmt );


    return 0;
    }


    // The member function getDbData(...) is a member
    // function for extracting information from the database.
    //
    // Returns the data from a table using a query of the form
    //
    // SELECT colName[0], colName[1],  ... FROM tableName WHERE colSpec[0] = colSpecVal[0] AND colSpec[1] = colSpecVal[1] ... ORDER BY orderBy
    //
    // as an array of arrays of strings. The use of strings to return all data is to reduce the complexity
    // of creating return arguments of variable type.
    //
    // The conversion of any returned data to integer or floating point values  is accomplished
    // using atoi or atof.
    //
    // The types are specified in colSpecType as strings, and one of INTEGER, REAL, or TEXT
    //
    // The orderBy string specifies the column to use for output ordering. If empty,
    // then that portion of the SQL statement is dropped.
    //
    // If distinctFlag is set, then the DISTINCT keyword is added to the SELECT statement
    //

	int getDbData(const string& tableName, const vector <string > & colName, const vector <string > & colSpec,
    const vector <string > & colSpecVal,  const vector <string > & colSpecType,  const string& orderBy,
    vector < vector< string> >& colData, bool distinctFlag = false)
    {
        int rc             = 0;
        sqlite3_stmt* stmt = 0;
	    ostringstream tmpStream;

        if(distinctFlag)
        {
		tmpStream << "SELECT DISTINCT ";
		}
		else
		{
		tmpStream << "SELECT ";
		}

        long colCount      = colName.size();
	    long colSpecCount  = colSpec.size();

		for(long i = 0; i < colCount; i++)
		{
		tmpStream << " " << colName[i].c_str() << " ";
		if(i < colCount-1) { tmpStream << ","; }
		}

        if(colSpecCount > 0)
        {
		tmpStream <<  "FROM " << tableName.c_str() << " WHERE ";
		}
		else
		{
		tmpStream <<  "FROM " << tableName.c_str();
		}

		for(long i = 0; i < colSpecCount; i++)
		{
		if(colSpecType[i].compare("TEXT") == 0)
	    {
		tmpStream << colSpec[i].c_str() << " = \'" << colSpecVal[i] << "\' ";
		if(i < colSpecCount-1) {tmpStream << " AND ";}
		}
		else
		{
		tmpStream << colSpec[i].c_str() << " = " << colSpecVal[i];
		if(i < colSpecCount-1) {tmpStream << " AND ";}
		}

		}

        if(orderBy.size() > 0)
        {
        tmpStream << " ORDER BY " << orderBy.c_str() << ";";
        }
        else
        {
        tmpStream << ";";
        }

        if(verboseDBflag) {printf("%s\n",tmpStream.str().c_str());}
        sqlUtil.createUTF8characterEncoding(tmpStream.str(),utf8chars);
        rc = sqlite3_prepare_v2(db, &utf8chars[0], -1, &stmt, 0);

        if(rc != SQLITE_OK)
        {
	 		printf(" SQ_ERROR           : %s \n",sqlite3_errmsg(db));
	 	    sqlite3_finalize( stmt );
            errorHandling();
	 		return 1;
        }

        // All values are read using the text function so that no matter how SQLite decides
        // to store the values they are converted to a common type, namely text.
        // The sqlite3_column_text() routine automatically converts the entries appropriately.
        //
        colData.clear();
        colData.resize(colCount);

        const char* sPtr;
        int    dbColType;
        string     realString;
        string     realStringTrimmed;
        string     nullString;

        while( sqlite3_step( stmt ) == SQLITE_ROW )
        {
        for(long k = 0; k < colCount; k++)
        {
        dbColType = sqlite3_column_type(stmt,k);
        switch (dbColType)
        {
        case SQLITE_INTEGER :
          sPtr    = (const char*)sqlite3_column_text(stmt, k);
          if(sPtr != 0){colData[k].push_back(sPtr);} else {colData[k].push_back(nullString);} break;
        case SQLITE_FLOAT   :
          sPtr    = (const char*)sqlite3_column_text(stmt, k);
          if(sPtr != 0)
          {realString.assign(sPtr); sqlUtil.removeScientificFormatTrailingZeros(realString,realStringTrimmed); colData[k].push_back(realStringTrimmed);}
          else {colData[k].push_back(nullString);} break;
        case SQLITE_TEXT    :
          sPtr    = (const char*)sqlite3_column_text(stmt, k);
          if(sPtr != 0){colData[k].push_back(sPtr);} else {colData[k].push_back(nullString);} break;
        default:
           colData[k].push_back(nullString);
        }
        }

        }
        sqlite3_finalize( stmt );
        return 0;
    }

//
//  This routine queries the database and captures the database
//  schema information into dbSchemaData, an instance of
//  map< string, map < string , string > >
//
//  This information is used for input/output checking when
//  extracting or updating information in the table.
//
//  dbSchemaData[tableName] returns a map< string, string>
//  whose index is the column name and the value is the column
//  data type (one of TEXT, REAL, or INTEGER).
//
//
    void updateDbSchemaData()
    {
    vector < string > tableNames;
    getDBtableNames(tableNames);

    vector < string > columnNames;
    vector < string > columnDataTypes;

    dbSchemaData.clear();

    map <string, string > tableSchema;
    for(vector<string>::iterator name = tableNames.begin(); name != tableNames.end(); ++name)
    {
    columnNames.clear();
    columnDataTypes.clear();
    tableSchema.clear();

    getDBcolumnInfo((*name),columnNames,columnDataTypes);

    tableSchema.clear();
    for(int i = 0; i < (int)columnNames.size(); i++)
    {
    tableSchema.insert(std::pair<string,string>(columnNames[i],columnDataTypes[i]));
    }
    dbSchemaData.insert(std::pair<string,map <string, string > > (*name, tableSchema));
    }
    }

    void printDbSchemaData()
    {
    std::map < string, map < string, string > >::iterator it;
    std::map < string, string >                     colData;
    string tableName;
    string colName;
    string colDataType;
    for(it =  dbSchemaData.begin(); it != dbSchemaData.end(); ++it)
    {
    	tableName = it->first;
    	printf("Table Name : %s \n",tableName.c_str());
    	colData   = it->second;
    	for(std::map < string , string >::iterator col = colData.begin(); col != colData.end(); ++col)
    	{
        colName     = col->first;
        colDataType = col->second;
        printf("%s  %s \n",colName.c_str(),colDataType.c_str());
    	}
    }
    }

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

    SimDb_Utilities sqlUtil;

	bool  connectionFlag;
	bool  existingDBflag;

	bool haltOnErrorFlag;
	bool defaultInitFlag;
	string  dataBaseName;

    vector < char > utf8chars;
    bool        verboseDBflag;
    bool          verboseFlag;

    sqlite3*               db;

    std::map< string, map  < string, string > > dbSchemaData;
};
#endif
