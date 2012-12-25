def _convert_mysqlquery_sqlite(sqlpath):
	''' Reads in the mysql .sql file and returns a sqlite compatible query

	'''
	noncompatible = ["ENGINE=MyISAM DEFAULT CHARSET=latin1", "SET character_set_client = @saved_cs_client;", "SET @saved_cs_client     = @@character_set_client;","SET character_set_client = utf8;","unsigned","auto_increment","MAX_ROWS=4294967295 AVG_ROW_LENGTH=40", " MAX_ROWS=4294967295","AVG_ROW_LENGTH=200","ENGINE=MyISAM AUTO_INCREMENT=46974 DEFAULT CHARSET=latin1"]
	comment_delim = ["--", "/*","*"]
		


	datatypes = {" char(" : " text ", " varchar(" : " text ", " tinytext " : " text ", " mediumtext " : " text ", " longtext " : " text ",
			" mediumblob " : " blob ", " longblob " : " blob ",
			" tinyint(" : " integer ", " int("  : " integer ", " smallint("  : " integer ", " mediumint("  : " integer ", " bigint("  : " integer ",
			" float " : " real ", " double("  : " real ", " decimal("  : " real ",
			" enum(" : " text ", " set(" : " text "}


	sqlquery = ""
	with open(sqlpath,"r") as sql:
		while True:
			line = sql.readline()
			if len(line) == 0:
				break
			for noncom in noncompatible:	
				line = line.replace(noncom,"")

			# this line is needed to remove the 'key' if it exists
			if "KEY" in line:
				line = ""
				
			# remove comment lines
			for comline in comment_delim:
				if comline in line:
					line = ""

			# remove mysql data types and replaces them with sqlite datatypes
			for dt in datatypes.keys():
				if dt in line:
					if "(" in dt:
						line = line.replace(dt,"*")
						dtindex = line.find("*")
						line = line[:dtindex] + datatypes[dt] + line[dtindex + len(datatypes[dt]):]
					else:
						line = line.replace(dt,datatypes[dt])


			sqlquery += line
				
	
	sqlquery = sqlquery.replace("\n","").replace("\r","")
	sqlquery = sqlquery.replace(",)",")")
	ind_last_comma = sqlquery.rfind(",")
	ind_last_par = sqlquery.rfind(")")

	if (ind_last_comma == -1) and (ind_last_par == -1) and ind_last_comma > ind_last_par:
		sqlquery = sqlquery[:ind_last_comma] + sqlquery[ind_last_comma + 1:]
	return sqlquery	


