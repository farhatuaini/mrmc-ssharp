#!/bin/bash -u

#Deleting files with the following patterns
DELETE_FILE_PATTERNS=".*~ .*\.bak .*\.backup .*\.out .*\.nav .*\.toc .*\.aux .*\.snm .*\.bbl .*\.blg .*\.log .*\.dvi .*\.ps";

#Internal variables and functions
clean_line() {
	#Erase the "${1}/" string from the output
	MSG_LENGTH=${#1};
	BSP_STRING="";
    BLANK_STRING="";
	let "NUMBER_OF_BSP = ${MSG_LENGTH} + 1";
	for (( i=1 ; i <= ${NUMBER_OF_BSP} ; i++ )) ; do
		BSP_STRING="${BSP_STRING}\b";
		BLANK_STRING="${BLANK_STRING} ";
	done
	echo -n -e "${BSP_STRING}${BLANK_STRING}${BSP_STRING}";
}

#Get the list of all files here
MSG="Get the list of all files";
echo -n ${MSG};
ALL_FILES=.tmp_files;
rm -f ${ALL_FILES};
find . > ${ALL_FILES};
clean_line "${MSG}";

#Select the files we want to delete
FILES_TO_DELETE=.tmp_delete_files;
rm -f ${FILES_TO_DELETE};
for pattern in ${DELETE_FILE_PATTERNS}
do
	MSG="Processing pattern: "${pattern};
	echo -n ${MSG};
	cat ${ALL_FILES} | grep "${pattern}$" >> ${FILES_TO_DELETE};
	clean_line "${MSG}";
done

#Removing the tmp files
MSG="Removing the tmp files";
echo -n ${MSG};
DELETE_FILES_LIST=`eval cat ${FILES_TO_DELETE}`
for one_file in ${DELETE_FILES_LIST}; do
	rm -f ${one_file}
done

rm -f ${ALL_FILES} ${FILES_TO_DELETE};
clean_line "${MSG}";

