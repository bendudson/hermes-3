#!/usr/bin/env bash

function doit() {
    dir=$1
    
    for d in $(git ls-tree --full-name HEAD --name-only  -r -- $dir|grep -o '.*/'|sort -u) ; do
	doit_f $d
    done
}

function doit_f() {
    dir=$1
    dirn=${dir//\//_}

    echo "list(APPEND copy_dir_${dirn}files"
    #git ls-tree --full-name HEAD --name-only  -r -- $dir
    for fn in $(git ls-tree --full-name HEAD --name-only -- $dir) ; do
	test -f $fn || continue
	echo "    $fn"
    done
    echo "    )"

    cat <<EOF

foreach(fn \${copy_dir_${dirn}files})
   list(APPEND copy_dir_${dirn}files_abs \${CMAKE_SOURCE_DIR}/\${fn})
endforeach()
message(INFO \${copy_dir_${dirn}files_abs})

add_custom_command(TARGET \${PROJECT_NAME} POST_BUILD
                   COMMAND \${CMAKE_COMMAND} -E make_directory
		   \$<TARGET_FILE_DIR:\${PROJECT_NAME}>/$dir)

add_custom_command(TARGET \${PROJECT_NAME} POST_BUILD
                   COMMAND \${CMAKE_COMMAND} -E copy_if_different
                   \${copy_dir_${dirn}files_abs} \$<TARGET_FILE_DIR:\${PROJECT_NAME}>/$dir)
EOF
}

outfile=${0/.bash/}

doit examples > $outfile

doit json_database >> $outfile

doit tests >> $outfile
