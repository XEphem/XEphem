
This application generates lists of eclipses via Inex and Saros series.

When modifying libastro/eclipse.c code, this application can be used to detect changes in calculated eclipse date/time or eclipse path.

Eclipse lists are saved in the data/ subdirectory.
The generated lists are grouped by Saros series or sorted by eclipse date.
The makefile default is generating Saros series 0 thru 180.

See https://eclipse.gsfc.nasa.gov/SEsaros/SEsaros.html for comparison data.

Examples:

# Build the application.
make

# Build the application and generate data.
make data

# Build the application, generate data, and run tests.
make test

# Update the test data hashes.
make hash

# Remove the built application.
make clean

# Remove the created test data.
make cleantest
