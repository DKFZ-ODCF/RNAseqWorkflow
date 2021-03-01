#!/usr/bin/env bash

IFS=$'\n\t'

set -vx

function echo_run (){
        local COMMAND="$*"
        echo $COMMAND
        if [[ "$TEST_RUN" != true ]]
        then
                eval $COMMAND
                if [[ $? -ne 0 ]]
                then
                  echo "Error running '$COMMAND'" && exit
                fi
        fi
}

function remove_directory () {
        local dir_to_remove="$1"
        if [ -d "$dir_to_remove" ]
        then
                echo "#DIRECTORY: \"$dir_to_remove\" exists... removing" 1>&2
		        rm -r "$dir_to_remove"
        else
                echo "#DIRECTORY: \"$dir_to_remove\" does not exist (or isn't a directory)... skipping" 1>&2
        fi
}

function remove_file () {
        local file_to_remove="$1"
        if [[ -e "$file_to_remove" && ! -d "$file_to_remove" ]]
        then
                echo "#FILE: \"$file_to_remove\" exists... removing" 1>&2
                rm "$file_to_remove"
        else
                echo "#FILE: \"$file_to_remove\" does not exist (or isn't a file)... skipping" 1>&2
        fi
}

# make directories if they do not exist
function make_directory () {
        local dir_to_make="$1"
        if [ -d "$dir_to_make" ]
        then
                echo "#DIRECTORY: \"$dir_to_make\" exists... skipping" 1>&2
        else
		echo "#DIRECTORY: \"$dir_to_make\" does not exist... creating" 1>&2
                mkdir -p "$dir_to_make"
        fi
}

# file check deaths
function check_or_die (){
        local file_to_check="$1"
        local what_is_running="$2"
        PWD=eval pwd
	    if [[ -f "$PWD/$file_to_check" ]]
            then
		    echo "#FOUND FILE: \"$PWD/$file_to_check\" for process \"$what_is_running\"... continuing run mode" 1>&2
	    elif [[ "$TEST_RUN" == true ]]
	    then
		    echo "#FILE NOT CHECKED: \"$PWD/$file_to_check\" for process \"$what_is_running\"... continuing test mode" 1>&2
	    else
		    echo "#ERROR: file not found: \"$PWD/$file_to_check\" for process \"$what_is_running\"... exitting!" 1>&2
		    exit 1
	    fi
}

function check_executable (){
        local exe_to_check="$1"
        local exe_with_path="$(which $exe_to_check)"
	    if [[ -x "$exe_with_path" ]]
            then
		    echo "#EXE FILE: \"$exe_to_check\" with path \"$exe_with_path\"... continuing run mode" 1>&2
	    elif [[ "$TEST_RUN" == true ]]
	    then
		    echo "#EXE FILE NOT CHECKED: \"$exe_to_check\" with path \"$exe_with_path\"... continuing test mode" 1>&2
	    else
		    echo "#ERROR: exe file not found (or executable): \"$$exe_to_check\" with path \"$exe_with_path\"... exitting!" 1>&2
		    exit 1
	    fi
}

function check_text_or_die (){
        local file_to_check="$1"
        local text_to_check="$2"
	local grep_text=`grep ${text_to_check} ${file_to_check}`
            if [[ ${#grep_text} -ge 5 ]]
            then
                    echo "#CHECK TEXT OR DIE: \"$file_to_check\" contains text \"$text_to_check\"... continuing run mode" 1>&2
	    elif [[ "$TEST_RUN" == true ]]
            then
                    echo "#TEXT NOT CHECKED: \"$file_to_check\" contains text \"$text_to_check\"... continuing test mode" 1>&2            
	    else
                    echo "#ERROR: \"$file_to_check\" DOES NOT contains text \"$text_to_check\"... exitting!" 1>&2
                    exit 1
            fi
}

function check_text_and_die (){
        local file_to_check="$1"
        local text_to_check="$2"
	local remedy="$3"
        local grep_text=`grep ${text_to_check} ${file_to_check}`
            if [[ ${#grep_text} -ge 5 ]]
            then
                    echo "#ERROR: \"$file_to_check\" contains text \"$text_to_check\"... $remedy ... exitting!" 1>&2
                    exit 1
            elif [[ "$TEST_RUN" == true ]]
            then
                    echo "#TEXT NOT CHECKED: \"$file_to_check\" contains text \"$text_to_check\"... continuing test mode" 1>&2
            else
                    echo "#CHECK TEXT AND DIE: \"$file_to_check\" DOES NOT contains text \"$text_to_check\"... continuing run mode" 1>&2
            fi
}

