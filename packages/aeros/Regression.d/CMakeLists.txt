# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
#add_library (Hello hello.cxx) 
add_custom_command(
    OUTPUT copy_files.pl
    COMMAND ${REG_SRC_DIR}/copy_files.pl $REG_SRC_DIR $EXEC_DIR
#   ${CMAKE_CURRENT_SOURCE_DIR}/schema.txt.in
#   ${CMAKE_CURRENT_SOURCE_DIR}/generated.h
    DEPENDS copy_files.pl
    COMMENT "copying files..."
)

