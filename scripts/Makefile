F03=gfortran
F03FLAGS=-O3 -g -fopenmp -Wall -frecord-marker=4 -ffree-line-length-none
DEBUGFLAGS = -g -fbacktrace -fbounds-check

FFLAGS = $(F03FLAGS) $(DEBUGFLAGS)

FFTWLIB= -I/$(FFTW_ROOT)/include -lfftw3 -lfftw3_omp




default: all

all: get_clump_ID_of_orphans \
		count_duplicate_orphans \
		get_halostats \
		get_correlation \
		get_smf \
		get_orphan_image_data

########################
# Execs
########################

get_clump_ID_of_orphans: quick_sort.o get_clump_ID_of_orphans.o io_module.o constants_and_parameters.o
	$(F03) $(FFLAGS) $^ -o $@

count_duplicate_orphans: quick_sort.o count_duplicate_orphans.o io_module.o constants_and_parameters.o
	$(F03) $(FFLAGS) $^ -o $@

get_halostats: quick_sort.o get_halostats.o io_module.o density_projection.o
	$(F03) $(FFLAGS) $^ -o $@

get_correlation: get_correlation.o io_module.o constants_and_parameters.o SMHM_relations.o
	$(F03) $(FFLAGS) $(FFTWLIB) $^ -o $@

get_smf: get_smf.o io_module.o constants_and_parameters.o SMHM_relations.o density_projection.o
	$(F03) $(FFLAGS) $^ -o $@

get_orphan_image_data: get_orphan_image_data.o io_module.o constants_and_parameters.o density_projection.o
	$(F03) $(FFLAGS) $^ -o $@



########################
# Modules
########################

constants_and_parameters.o: constants_and_parameters.f03
	$(F03) $(FFLAGS) -c constants_and_parameters.f03

quick_sort.o: constants_and_parameters.o quick_sort.f03
	$(F03) $(FFLAGS) -c quick_sort.f03

io_module.o: constants_and_parameters.o io_module.f03
	$(F03) $(FFLAGS) -c io_module.f03

SMHM_relations.o: constants_and_parameters.o SMHM_relations.f03
	$(F03) $(FFLAGS) -c SMHM_relations.f03

density_projection.o: constants_and_parameters.o io_module.o density_projection.f03
	$(F03) $(FFLAGS) -c density_projection.f03


########################
# Main program files
########################

count_duplicate_orphans.o: count_duplicate_orphans.f03 io_module.o quick_sort.o constants_and_parameters.o
	$(F03) $(FFLAGS) -c count_duplicate_orphans.f03

get_clump_ID_of_orphans.o: get_clump_ID_of_orphans.f03 io_module.o quick_sort.o constants_and_parameters.o
	$(F03) $(FFLAGS) -c get_clump_ID_of_orphans.f03

get_halostats.o: io_module.o quick_sort.o density_projection.o constants_and_parameters.o get_halostats.f03
	$(F03) $(FFLAGS) -c get_halostats.f03

get_smf.o: get_smf.f03 SMHM_relations.o io_module.o constants_and_parameters.o SMHM_relations.o density_projection.o
	$(F03) $(FFLAGS) -c get_smf.f03

get_correlation.o: get_correlation.f03 io_module.o constants_and_parameters.o SMHM_relations.o
	$(F03) $(FFLAGS) $(FFTWLIB) -c get_correlation.f03

get_orphan_image_data.o: get_orphan_image_data.f03 constants_and_parameters.o io_module.o
	$(F03) $(FFLAGS) -c get_orphan_image_data.f03


clean:
	rm -f *.o *.mod
