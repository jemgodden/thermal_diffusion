/*
 Title:   Solving the Thermal Diffusion Equation
 Author:  Jeremy Godden
*/

static const char * VERSION  = "1.0.1";
static const char * REV_DATE = "12-Mar-2020";

/*
 This program  has the ability to model how the temperature of a medium changes as a nuclear waste rod diffuses heat
 into its surroundings, making use of the cylindrical symmetry of the rod. This is predominantly an attempt to recreate
 Dr Olsen-Kettle's graph, referenced in the accompanying report, using her method.
 This program can also model how temperature changes radially in an egg as it is surrounded, and cooked, by boiling
 water, making use of the spherical symmetry, assuming it is spherical, of the egg. This uses the same method as
 used for modelling a nuclear waste rod, but is adapted for this scenario.

 An example of the command line arguments would be:
    ./a.out --egg --mass 60.0 --temperature 277.0 --file 60g_egg_data.txt
 This will model an egg of mass 60g being boiled, from an initial temperature of 277K. The data will be printed to the
 file 60g_egg_data.txt.
 Checks will be made on all command line argument values to make sure they are suitable.
 More information on these command line arguments can be found in the help() function below.

 All command line arguments are optional, apart from which model you would like to do. Default values will be set if
 none are specified. For the nuclear waste rod, a spatial resolution of 1%, and a time step of 0.1 years will be given.
 The output file will be rod_data.txt. For the egg model, it will be soft boiled, a spatial resolution of 1%, a time
 step of 1 second, a mass of 57g and an initial temperature of 293k (room temperature) will be given. The output file
 will be egg_data.txt.

 No input file is required as all parameters of the model can be specified via the command line arguments.

 When modelling a nuclear waste rod, the temperature at each radius will be printed to the output file at pre-determined
 times, in order to recreate Dr Olsen-Kettle's graph. When modelling an egg, the temperature at each radius will be
 printed to the output file at minute intervals.
 These files will then used to plot the temperature changes during the simulation, through GNUplot. GNUplot is accessed
 remotely through the terminal in this program.
*/

/* Packages used throughout the program. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

/* Defined values shared in both models. */
#define MATRIX_ELEMENT ((j * mesh_size) + i) /* Finds a specific element in a matrix. */
#define R_I ((i + 1) * delta_r) /* cm. Used to find the distance along the mesh when looping through i. */
#define ONE_HUNDRED_PERCENT 100.0
#define MINUTE 60.0 /* s */
#define ONE_THIRD (1.0 / 3.0) /* Double value of 1/3 used in powers. */
#define FILE_NAME_BUFFER 50 /* Buffer for length of file name being used. */
#define GNUPLOT "gnuplot -persist" /* Used to pipe to GNUplot to plot data. */

/* Defined values used in modelling a nuclear waste rod. */
#define TEMP_ENV 300.0 /* K */
#define TEMP_ROD 1.0 /* K */
#define TAU 100.0 /* s^-1 */
#define A 25.0 /* cm */
#define MAX_R 100.0 /* cm */
#define THERMAL_DIFFUSIVITY 2.0e7 /* cm^2 yr^-1 */
#define TIME_K (k * rod->time_step) /* yrs */
#define ROD_NO_OF_PLOTS 5 /* Number of plots used when modelling nuclear waste rod. */
#define DEFAULT_ROD_TIME_STEP 0.1 /* yrs */ /* Gives 1000 total steps, as done by Dr Olsen-Kettle. */
#define DEFAULT_ROD_RESOLUTION 1.010101 /* % */ /* Gives a mesh size of 99, as done by Dr Olsen-Kettle. */
#define DEFAULT_ROD_OUTPUT_FILE "rod_data.txt"

/* Defined values used in modelling an egg being boiled. */
#define ROOM_TEMP 294.0 /* K */
#define TEMP_WATER 373.0 /* K */
#define TEMP_COOKED 336.0 /* K */
#define YOLK_PERCENTAGE_MASS 0.31
#define WHITE_PERCENTAGE_MASS 0.58
#define SHELL_PERCENTAGE_MASS 0.11
#define YOLK_HEAT_CAPACITY 2.7 /* J g^-1 K^-1 */
#define WHITE_HEAT_CAPACITY 3.7 /* J g^-1 K^-1 */
#define SHELL_HEAT_CAPACITY 0.89 /* J g^-1 K^-1 */
#define YOLK_CONDUCTIVITY 3.4e-3 /* W cm^-1 K^-1 */
#define WHITE_CONDUCTIVITY 5.4e-3 /* W cm^-1 K^-1 */
#define SHELL_CONDUCTIVITY 2.3e-2 /* W cm^-1 K^-1 */
#define YOLK_DENSITY 1.032 /* g cm^-3 */
#define WHITE_DENSITY 1.038 /* g cm^-3 */
#define SHELL_DENSITY 2.3 /* g cm^-3 */
#define EGG_MAX_PLOTS 50 /* Maximum number of plots for modelling an egg. */
#define DEFAULT_EGG_TIME_STEP 0.1 /* s */
#define DEFAULT_EGG_RESOLUTION 0.1 /* % */
#define DEFAULT_EGG_MASS 57.0 /* g */
#define DEFAULT_EGG_OUTPUT_FILE "egg_data.txt"

/* Defined value constraints used in allowing and disallowing given arguments. */
#define MAX_RESOLUTION 10.0 /* % */
#define MIN_RESOLUTION 0.0 /* % */
#define MAX_ROD_TIME_STEP 1.0 /* yrs */
#define MAX_EGG_TIME_STEP 10.0 /* s */
#define MIN_TIME_STEP 0.0 /* s */
#define MAX_EGG_MASS 200.0 /* g */
#define MIN_EGG_MASS 1.0 /* g */
#define MIN_EGG_TEMP 0.0 /* K */

/* Flag set by ‘--verbose’. */
static int verbose_flag;

/* All error codes used within the program. */
typedef enum error_code {
    NO_ERROR = 0,
    HELP_CALLED = 1,
    INVALID_ARGS = 2,
    MEMORY_ERROR = 3,
    FILE_OPEN_ERROR = 4,
    PIPE_ERROR = 5,
    GSL_USE_ERROR = 6
} Error_code;

typedef int Error;


/* Information on the matrix being solved for when modelling the nuclear waste rod. */
typedef struct matrix {
    /* Pointers to gsl library vectors that are alloced and used in solving equation. */
    gsl_vector *super_diagonal;
    gsl_vector *diagonal;
    gsl_vector *sub_diagonal;
} Matrix;

/* All information used when modelling a nuclear waste rod. */
typedef struct rod {
    char file_name[FILE_NAME_BUFFER];
    double resolution;
    size_t mesh_size;
    double time_step;
    double delta_r;
    double s;
    /* Pointers to memory allocated structure. */
    Matrix *matrix;
    /* Pointers to memory allocated gsl library vectors used in the matrix equation solver. */
    gsl_vector *temp_vector;
    gsl_vector *rhs_vector;
} Rod;

/* All information used when modelling an egg being boiled. */
typedef struct egg {
    char file_name[FILE_NAME_BUFFER];
    /* Option for whether the egg is soft boiled or hard boiled.  */
    bool hard_boil;
    double mass;
    double resolution;
    double time_step;
    double initial_temperature;
    double yolk_radius;
    double white_radius;
    double radius;
    size_t mesh_size;
    size_t cooking_point;
    /* Values calculated once and saved in structure for future usage. */
    double delta_r;
    size_t yolk_white_boundary;
    size_t white_shell_boundary;
    double yolk_s;
    double white_s;
    double shell_s;
    /* Pointers to memory allocated structure. */
    Matrix *matrix;
    /* Pointer to memory allocated gsl library vectors used in the matrix equation solver. */
    gsl_vector *temp_vector;
    gsl_vector *rhs_vector;
} Egg;


/* Print help about command line arguments to the user. */
void print_help() {
    fprintf(stderr, "\nPlease enter in command lines like this:\n"
            "./a.out -OPTION [-s TIME_STEP] [-r RESOLUTION] [-m EGG_MASS] [-t EGG_TEMPERATURE] [-f OUTPUT_FILE]\n");
    fprintf(stderr, "All arguments in brackets are optional and will be set to default values if no value is given.\n");
    fprintf(stderr, "An example would be:\n"
            "\t./a.out -e -s 0.5 -m 60.0 -f 60g_egg_data.txt\n");
    fprintf(stderr, "Options, which must be chosen first:\n"
            "\t-h: Shows help with command line arguments.\n"
            "\t-n: Models a nuclear waste rod.\n"
            "\t-e: Models and egg being boiled.\n");
    fprintf(stderr, "Optional arguments:\n"
            "\t-s TIME_STEP, is the step between each new position calculation, in seconds.\n"
            "\t-r RESOLUTION, is the spatial resolution of the model, in percentage.\n"
            "\t-r HARD_BOIL, is the option to hard boil the egg, instead of soft boiling the egg, by default. This is only applicable when modelling an egg.\n"
            "\t-m EGG_MASS, is the mass of the egg being modelled, in grams. This is only applicable when modelling an egg.\n"
            "\t-t EGG_TEMPERATURE, is the initial temperature of the egg being modelled, in kelvin. This is only applicable when modelling an egg.\n"
            "\t-f OUTPUT_FILE, is the name of the file that the data for the simulation will be printed to.\n\n");
}


/* Print an error message due to invalid command line arguments. */
Error invalid_args_print(const char *message, const char *option, const char *value) {
    fprintf(stderr, "\nInvalid command line arguments were given. %s\n", message);
    fprintf(stderr, "There error came in option: %s, and may be due to the value given: %s\n", option, value);
    print_help();
    return INVALID_ARGS;
}

/* Print an error message due to inability to allocate memory. */
Error memory_error_print() {
    fprintf(stderr, "\nMemory could not be allocated.\n");
    return MEMORY_ERROR;
}

/* Print an error message due to inability to pipe to GNUplot. */
Error pipe_error_print() {
    fprintf(stderr, "\nPipe to GNUplot could not be opened.\n");
    return PIPE_ERROR;
}

/* Print an error message due to inability to open a file. */
Error file_open_error_print(const char *file_name) {
    fprintf(stderr, "\nFile %s could not be opened.\n", file_name);
    return FILE_OPEN_ERROR;
}

/* Print an error message due a gsl function returning an error. */
Error gsl_error_print(const char *message, const int error) {
    fprintf(stderr, "\nThere was an error using the GSL library. %s\n", message);
    fprintf(stderr, "The driver returned error code %d.\n", error); /* Outlines error code from gsl library for the user. */
    return GSL_USE_ERROR;
}


/* Safely convert a double into an int, using ceil function to round the value up to the nearest integer. */
int ceil_double_to_int(double double_value) {
    double_value = ceil(double_value);
    /* Convert to an int by casting. */
    return (int)double_value;
}

/* Safely convert a double into a size_t, using floor function to round the value up to the nearest integer. */
size_t floor_double_to_size_t(double double_value) {
    double_value = floor(double_value);
    /* Convert to a size_t by casting. */
    return (size_t)double_value;
}

/* Convert command line argument string to a double. */
Error get_arg_double(const char *arg_value, double *value, const char *option) {
    char *end_ptr;
    *value = strtod(arg_value, &end_ptr);

    /* Checks that there are no extra characters after an argument given. */
    if (*end_ptr != '\0') {
        return invalid_args_print("Invalid number given.", option, arg_value);
    }

    return NO_ERROR;
}

/* Check if gsl_vector was successfully allocated and safely frees if so. */
void safe_gsl_vector_free(gsl_vector *vector) {
    if (vector != NULL) {
        gsl_vector_free(vector);
        vector = NULL; /* For safety. */
    }
}

/* Check if a pointer to a variable was successfully allocated and safely frees if so. */
void safe_free(void *ptr) {
    if (ptr != NULL) {
        free(ptr);
        ptr = NULL; /* For safety. */
    }
}


/* Free all memory to do with the matrix structure. */
void free_matrix(Matrix *matrix) {
    safe_gsl_vector_free(matrix->super_diagonal);
    safe_gsl_vector_free(matrix->diagonal);
    safe_gsl_vector_free(matrix->sub_diagonal);
    safe_free(matrix);
}

/* Create the matrix structure. */
Error create_matrix(Matrix *matrix, const size_t mesh_size) {
    /* Allocating memory for vectors in structure. */
    matrix->super_diagonal = gsl_vector_alloc(mesh_size - 1);
    matrix->diagonal = gsl_vector_alloc(mesh_size);
    matrix->sub_diagonal = gsl_vector_alloc(mesh_size - 1);
    if (matrix->super_diagonal == NULL || matrix->diagonal == NULL || matrix->sub_diagonal == NULL) {
        return memory_error_print();
    }

    return NO_ERROR;
}

/* Give each point in the mesh an initial temperature. */
void setup_initial_temp_vector(gsl_vector *temp_vector, const double initial_temp) {
    gsl_vector_set_all(temp_vector, initial_temp);
}

/* Save the temperature for all radii of the model into an array of saved arrays to be plotted later. */
void save_data(const gsl_vector *temp_vector, double *saved_temps, const size_t mesh_size, const int j) {
    for (size_t i=0; i<mesh_size; i++) {
        /* Store the saved data arrays in a matrix. */
        saved_temps[MATRIX_ELEMENT] = gsl_vector_get(temp_vector, i);
    }
}

/* Print all saved data, from the model, to a file. */
void file_print_data(FILE *f, const double *saved_temps, const size_t mesh_size, const double delta_r, const int no_plots) {
    for (size_t i=0; i<mesh_size; i++) {
        /* Print radius of each temperature first. */
        fprintf(f, "%lf", R_I);
        for (int j=0; j<no_plots; j++) {
            /* Print temperature at each radius for each saved array of data. */
            fprintf(f, "\t%lf", saved_temps[MATRIX_ELEMENT]);
        }
        fprintf(f, "\n");
    }
}

/* Solve the matrix equation for the model. */
Error solve_matrix_equation(Matrix *matrix, gsl_vector *rhs_vector, gsl_vector *temp_vector) {
    Error error;

    /* Using the gsl library in order to solve the matrix equation. */
    error = gsl_linalg_solve_tridiag(matrix->diagonal, matrix->super_diagonal, matrix->sub_diagonal, rhs_vector, temp_vector);
    if (error != NO_ERROR) {
        /* Returns error to user if gsl function returns an error. */
        return gsl_error_print("", error);
    }

    return NO_ERROR;
}


/* Free all memory to do with the rod structure. */
void free_rod(Rod *rod) {
    if (rod->matrix != NULL) {
        /* Check if pointer to matrix structure has been allocated, before attempting to free it. */
        free_matrix(rod->matrix);
    }
    safe_gsl_vector_free(rod->temp_vector);
    safe_gsl_vector_free(rod->rhs_vector);
    safe_free(rod);
}

/* Create the diagonals of the matrix being solved. */
void rod_setup_matrix_diags(Rod *rod) {
    /* Setting all values of the diagonal matrix to be the same. */
    gsl_vector_set_all(rod->matrix->diagonal, 1 + (2 * rod->s));
    /* First term of diagonal is different to rest, so set after all other terms. */
    gsl_vector_set(rod->matrix->diagonal, 0, 1 + rod->s + (rod->s / 2));

    /* Creates sub and super diagonals one value shorter than diagonal. */
    for (size_t i=0; i<rod->mesh_size - 1; i++) {
        gsl_vector_set(rod->matrix->super_diagonal, i, - rod->s - (rod->s / (2 * (i + 1))));
        gsl_vector_set(rod->matrix->sub_diagonal, i, - rod->s + (rod->s / (2 * (i + 2))));
    }
}

/* Create the rod structure and populate it with information. */
Error create_rod(Rod *rod) {
    Error error;
    /* Finds size of the mesh being used. */
    rod->mesh_size = floor_double_to_size_t(ONE_HUNDRED_PERCENT / rod->resolution);

    /* Allocating memory for vectors in rod structure. */
    rod->temp_vector = gsl_vector_alloc(rod->mesh_size);
    rod->rhs_vector = gsl_vector_alloc(rod->mesh_size);

    /* Allocating memory for the matrix structure. */
    rod->matrix = malloc(sizeof(Matrix));
    if (rod->temp_vector == NULL || rod->rhs_vector == NULL || rod->matrix == NULL) {
        return memory_error_print();
    }

    /* Find information used in modelling the rod and save it in the rod structure for later use. */
    rod->delta_r = MAX_R / rod->mesh_size;
    rod->s = (THERMAL_DIFFUSIVITY * rod->time_step) / (pow(rod->delta_r, 2));

    setup_initial_temp_vector(rod->temp_vector, TEMP_ENV);

    error = create_matrix(rod->matrix, rod->mesh_size);
    if (error != NO_ERROR) {
        return error;
    }

    /* Give matrix structure vectors correct values. */
    rod_setup_matrix_diags(rod);

    return NO_ERROR;
}

/* Find the temperature contributed to the system, at each radius, due to the rod's radiation. */
double rod_source_value(const Rod *rod, const int k, const double r) {
    if (r <= A) {
        /* Only returns a non-zero function when within the radius of the rod. */
        return (exp(- TIME_K / TAU) * (TEMP_ROD / pow(A, 2)));
    }
    else {
        return 0.0;
    }
}

/* Set up the right hand side of the matrix equation being solved. */
void rod_setup_rhs_equation(const Rod *rod, const int k) {
    /* Combines all relevant vectors and values into one vector, which represents the RHS of the matrix equation. */
    for (size_t i=0; i<rod->mesh_size; i++) {
        if (i == rod->mesh_size - 1) {
            /* Final point in mesh is set up differently due to boundary conditions. */
            gsl_vector_set(rod->rhs_vector, i, gsl_vector_get(rod->temp_vector, i) -
                    ((- rod->s - (rod->s / (2 * rod->mesh_size))) * TEMP_ENV));
        }
        else {
            gsl_vector_set(rod->rhs_vector, i, gsl_vector_get(rod->temp_vector, i) + (THERMAL_DIFFUSIVITY *
                    rod->time_step * rod_source_value(rod, k, (i + 1) * rod->delta_r)));
        }
    }
}

/* Print all saved data of rod model to a file, that will be plotted using GNUplot. */
Error rod_print_data(const Rod *rod, const double *saved_rod_temps) {
    FILE *f = fopen(rod->file_name, "w+");
    if (f==NULL){
        return file_open_error_print(rod->file_name);
    }

    /* Print background information at top of file. */
    fprintf(f, "# Version = %s, Revision date = %s\n", VERSION, REV_DATE);
    fprintf(f, "# Temperatures radially out from the centre of a nuclear waste rod, and it's surroundings, "
            "at different times.\n");
    fprintf(f, "# Radius(cm)\tInitial Temp(K)\tTemp(K)[1 yr]\tTemp(K)[10 yrs]\tTemp(K)[50 yrs]\tTemp(K)[100 yrs]\n");

    file_print_data(f, saved_rod_temps, rod->mesh_size, rod->delta_r, ROD_NO_OF_PLOTS);

    fclose(f);
    return NO_ERROR;
}

/* Plot data obtained when simulating a nuclear waste rod, using GNUplot. */
Error rod_plot(const Rod *rod, const int *years) {
    FILE *gp;

    gp = popen(GNUPLOT, "w");
    if (gp==NULL) {
        return pipe_error_print();
    }

    const char *prefix = "plot ";

    fprintf(gp, "set title \"Temperature outwards from centre of a nuclear waste rod\" font \"Arial, 20\"\n");
    fprintf(gp, "set xlabel \"Radius (cm)\" font \"Arial, 22\"\n");
    fprintf(gp, "set ylabel \"Temperature (K)\" font \"Arial, 22\"\n");
    fprintf(gp, "set ylabel offset -2.5,0\n");
    fprintf(gp, "set xrange [ * : * ] noreverse writeback\n");
    fprintf(gp, "set yrange [ * : * ] noreverse writeback\n");
    fprintf(gp, "set tics font \"Arial, 15\"\n");
    fprintf(gp, "set key font \"Arial, 18\"\n");

    for (int i=0; i<ROD_NO_OF_PLOTS; i++) {
        /* Loop input to GNUplot in order to plot all data sets. */
        fprintf(gp, "%s\"%s\" using 1:%d title \'%d years\' with lines lw 4", prefix, rod->file_name, i + 2, years[i]);
        prefix = ", ";
    }
    fprintf(gp, "\n");

    pclose(gp);
    return NO_ERROR;
}

/* Find number of time steps required to get to specified years. */
void find_years_in_steps(const Rod *rod, const int *years, double *years_in_steps) {
    for (int i=0; i<ROD_NO_OF_PLOTS; i++){
        years_in_steps[i] = ceil_double_to_int(years[i] / rod->time_step);
    }
}

/* Model the nuclear waste rod. */
Error rod_model(Rod *rod, double *saved_rod_temps, const int *years) {
    Error error;

    /* Find number of steps at which data should be saved. */
    double years_in_steps[ROD_NO_OF_PLOTS];
    find_years_in_steps(rod, years, years_in_steps);

    /* Count to keep a track of number of sets of data saved, used to append data to correct point in saved array. */
    int saved_data_count = 0;
    /* Copy initial temperature to saved array. */
    save_data(rod->temp_vector, saved_rod_temps, rod->mesh_size, saved_data_count);
    saved_data_count++;

    for (int k=0; k<years_in_steps[4]; k++) {
        /* Loop over all time steps, finding the new temperatures each time. */
        rod_setup_rhs_equation(rod, k);
        error = solve_matrix_equation(rod->matrix, rod->rhs_vector, rod->temp_vector);
        if (error != NO_ERROR) {
            return error;
        }

        /* Save temperatures at specific times. */
        if (k == years_in_steps[1] || k == years_in_steps[2] || k == years_in_steps[3]) {
            /* Copy temperatures to array of saved values. */
            save_data(rod->temp_vector, saved_rod_temps, rod->mesh_size, saved_data_count);
            saved_data_count++;
        }
    }
    /* Copy final temperature to saved array. */
    save_data(rod->temp_vector, saved_rod_temps, rod->mesh_size, saved_data_count);

    /* Print data of rod temperatures at specific times to a file. */
    error = rod_print_data(rod, saved_rod_temps);
    if (error != NO_ERROR) {
        safe_free(saved_rod_temps);
        return error;
    }

    return NO_ERROR;
}

/* Set up information on rod and model it. */
Error rod_setup_model(Rod *rod) {
    Error error;

    error = create_rod(rod);
    if (error != NO_ERROR) {
        return error;
    }

    printf("A nuclear waste rod will be modelled.\n");

    /* Create array of values corresponding to year of each data set. */
    int years[ROD_NO_OF_PLOTS] = {0, 1, 10, 50, 100};

    /* Initialise and allocate memory for a large array to store values of temperature later used for plotting. */
    double *saved_rod_temps = malloc(sizeof(double) * rod->mesh_size * ROD_NO_OF_PLOTS);
    if (saved_rod_temps == NULL) {
        return memory_error_print();
    }

    /* Process of simulating the nuclear waste rod and recording data from it. */
    error = rod_model(rod, saved_rod_temps, years);
    if (error != NO_ERROR) {
        safe_free(saved_rod_temps);
        return error;
    }

    /* Free saved data as it is no longer needed. */
    safe_free(saved_rod_temps);

    /* Plot images of simulation for user to see. */
    error = rod_plot(rod, years);
    if (error != NO_ERROR) {
        return error;
    }

    printf("Data has been written to the file: %s.\n", rod->file_name);

    return NO_ERROR;
}

/* Set default values for modelling a rod. */
void give_rod_default_values(Rod *rod) {
    rod->time_step = DEFAULT_ROD_TIME_STEP;
    rod->resolution = DEFAULT_ROD_RESOLUTION;
    strcpy(rod->file_name, DEFAULT_ROD_OUTPUT_FILE);
}


/* Function to free all memory to do with the egg structure. */
void free_egg(Egg *egg) {
    if (egg->matrix != NULL) {
        /* Check if pointer to matrix structure has been allocated, before attempting to free it. */
        free_matrix(egg->matrix);
    }
    safe_gsl_vector_free(egg->temp_vector);
    safe_gsl_vector_free(egg->rhs_vector);
    safe_free(egg);
}

/* Create the diagonals of the matrix being solved. */
void egg_setup_matrix_diags(Egg *egg) {
    for (size_t i=0; i<egg->yolk_white_boundary; i++) {
        /* Give all matrix values corresponding to within the yolk radius the correct value. */
        gsl_vector_set(egg->matrix->diagonal, i, 1 + (2 * egg->yolk_s));
        gsl_vector_set(egg->matrix->super_diagonal, i, - egg->yolk_s - (egg->yolk_s / (i + 1)));
        if (i != 0) {
            /* Start the sub diagonal one term later than the diagonal in order to give it one less term overall. */
            gsl_vector_set(egg->matrix->sub_diagonal, i - 1, - egg->yolk_s + (egg->yolk_s / (i + 1)));
        }
    }

    for (size_t j=egg->yolk_white_boundary; j<egg->white_shell_boundary; j++) {
        /* Give all matrix values corresponding to within the yolk and white radii the correct value. */
        gsl_vector_set(egg->matrix->diagonal, j, 1 + (2 * egg->white_s));
        gsl_vector_set(egg->matrix->super_diagonal, j, -egg->white_s - (egg->white_s / (j + 1)));
        gsl_vector_set(egg->matrix->sub_diagonal, j - 1, - egg->white_s + (egg->white_s / (j + 1)));
    }

    for (size_t k=egg->white_shell_boundary; k<egg->mesh_size; k++) {
        /* Give all matrix values corresponding to within the white and shell radii the correct value. */
        gsl_vector_set(egg->matrix->diagonal, k, 1 + (2 * egg->shell_s));
        if (k != egg->mesh_size - 1) {
            /* End the super diagonal one term shorter than the diagonal. */
            gsl_vector_set(egg->matrix->super_diagonal, k, - egg->shell_s - (egg->shell_s / (k + 1)));
        }
        gsl_vector_set(egg->matrix->sub_diagonal, k - 1, - egg->shell_s + (egg->shell_s / (k + 1)));
    }
}

/* Find relevant information of the egg, to be used in the simulation. */
void find_egg_information(Egg *egg) {
    /* Find radius of all egg components. */
    egg->yolk_radius = pow(3 * YOLK_PERCENTAGE_MASS * egg->mass / (4 * M_PI * YOLK_DENSITY), ONE_THIRD);
    egg->white_radius = pow((pow(egg->yolk_radius, 3) + ((3 * WHITE_PERCENTAGE_MASS * egg->mass) / (4 * M_PI * WHITE_DENSITY))), ONE_THIRD);
    egg->radius = pow((pow(egg->white_radius, 3) + ((3 * SHELL_PERCENTAGE_MASS * egg->mass) / (4 * M_PI * SHELL_DENSITY))), ONE_THIRD);

    /* Find boundaries of egg components. */
    egg->yolk_white_boundary = floor_double_to_size_t((egg->yolk_radius / egg->radius) * egg->mesh_size);
    egg->white_shell_boundary = floor_double_to_size_t((egg->white_radius / egg->radius) * egg->mesh_size);

    egg->delta_r = egg->radius / egg->mesh_size;

    /* Calculate egg component S values, as specified in formulas, for use later in the program. */
    egg->yolk_s = (YOLK_CONDUCTIVITY * egg->time_step) / (YOLK_HEAT_CAPACITY * YOLK_DENSITY * pow(egg->delta_r, 2));
    egg->white_s = (WHITE_CONDUCTIVITY * egg->time_step) / (WHITE_HEAT_CAPACITY * WHITE_DENSITY * pow(egg->delta_r, 2));
    egg->shell_s = (SHELL_CONDUCTIVITY * egg->time_step) / (SHELL_HEAT_CAPACITY * SHELL_DENSITY * pow(egg->delta_r, 2));

    /* Check to see if egg should be hard boiled, and set appropriate cooking point of mesh. */
    if (egg->hard_boil == true) {
        egg->cooking_point = 0;
    }
    else {
        egg->cooking_point = egg->yolk_white_boundary;
    }
}

/* Create the egg structure and populate it with information. */
Error create_egg(Egg *egg) {
    Error error;
    /* Finds size fo mesh being used. */
    egg->mesh_size = floor_double_to_size_t(ONE_HUNDRED_PERCENT / egg->resolution);

    /* Allocating memory for vectors in egg structure. */
    egg->temp_vector = gsl_vector_alloc(egg->mesh_size);
    egg->rhs_vector = gsl_vector_alloc(egg->mesh_size);
    /* Allocating memory for the matrix structure. */
    egg->matrix = malloc(sizeof(Matrix));
    if (egg->temp_vector == NULL || egg->rhs_vector == NULL || egg->matrix == NULL) {
        return memory_error_print();
    }

    find_egg_information(egg);

    setup_initial_temp_vector(egg->temp_vector, egg->initial_temperature);

    /* Allocate correct sized memory for matrix structure vectors. */
    error = create_matrix(egg->matrix, egg->mesh_size);
    if (error != NO_ERROR) {
        return error;
    }

    egg_setup_matrix_diags(egg);

    return NO_ERROR;
}

/* Set up the right hand side of the matrix equation being solved for an egg model. */
void egg_setup_rhs_equation(Egg *egg) {
    /* Combines all relevant vectors and values into one vector, which represents the RHS of the matrix equation. */
    for (size_t i=0; i<egg->mesh_size; i++) {
        if (i == egg->mesh_size - 1) {
            /* Final point in mesh is set up differently due to boundary conditions. */
            gsl_vector_set(egg->rhs_vector, i, gsl_vector_get(egg->temp_vector, i) -
                    ((- egg->shell_s - (egg->shell_s / egg->mesh_size)) * TEMP_WATER));
        }
        else {
            gsl_vector_set(egg->rhs_vector, i, gsl_vector_get(egg->temp_vector, i));
        }
    }
}

/* Print all saved data of egg model to a file, that will be plotted using GNUplot. */
Error egg_print_data(const Egg *egg, const double *saved_egg_temps, const int saved_arrays_count, const double time) {
    FILE *f = fopen(egg->file_name, "w+");
    if (f==NULL){
        return file_open_error_print(egg->file_name);
    }

    /* Print background information at top of file. */
    fprintf(f, "# Version = %s, Revision date = %s\n", VERSION, REV_DATE);
    fprintf(f, "# Temperatures radially out from the centre of an egg, as it's being boiled, at different times.\n");
    /* Print mass of egg and time it took to cook. */
    fprintf(f, "# This egg had a mass of %.1lfg and it took %.1lf seconds to cook.\n", egg->mass, time);
    fprintf(f, "# Radius(cm)\tInitial Temp(K)");
    for (int k=0; k<(saved_arrays_count-1); k++) {
        fprintf(f, "\tTemp(K)[%d mins]", (k + 1));
    }
    fprintf(f, "\tFinal Temp(K)\n");

    file_print_data(f, saved_egg_temps, egg->mesh_size, egg->delta_r, saved_arrays_count);

    fclose(f);
    return NO_ERROR;
}

/* Plot data obtained when modelling an egg being boiled, using GNUplot. */
Error egg_plot(const Egg *egg, const int saved_arrays) {
    FILE *gp;

    gp = popen(GNUPLOT, "w");
    if (gp==NULL) {
        return pipe_error_print();
    }

    const char *prefix = "plot ";

    fprintf(gp, "set title \"Radial temperature of an egg in boiling water\" font \"Arial, 20\"\n");
    fprintf(gp, "set xlabel \"Radius (cm)\" font \"Arial, 22\"\n");
    fprintf(gp, "set ylabel \"Temperature (K)\" font \"Arial, 22\"\n");
    fprintf(gp, "set ylabel offset -2.5,0\n");
    fprintf(gp, "set xrange [ * : * ] noreverse writeback\n");
    fprintf(gp, "set yrange [ * : * ] noreverse writeback\n");
    fprintf(gp, "set tics font \"Arial, 15\"\n");
    fprintf(gp, "set key above\n");
    fprintf(gp, "set key font \"Arial, 15\"\n");

    for (int i=0; i<saved_arrays; i++) {
        /* Loop input to GNUplot in order to plot all data sets. */
        if (i == 0) {
            /* If first plot, label it with the start of the simulation. */
            fprintf(gp, "%s \"%s\" using 1:%d title \'Start\' with lines lw 4", prefix, egg->file_name, (i + 2));
            prefix = ", ";
        }
        else if (i == saved_arrays - 1) {
            /* If first plot, label it as the point at which the egg became cooked. */
            fprintf(gp, "%s \"%s\" using 1:%d title \'Cooked\' with lines lw 4", prefix, egg->file_name, (i + 2));
        }
        else {
            fprintf(gp, "%s \"%s\" using 1:%d title \'After %d minutes\' with lines lw 4", prefix, egg->file_name,
                    (i + 2), i);
            prefix = ", ";
        }
    }
    fprintf(gp, "\n");

    pclose(gp);
    return NO_ERROR;
}

/* Model the egg being boiled. */
Error egg_model(Egg *egg, double *saved_egg_temps, int *saved_data_count) {
    Error error;

    /* Start time at 0 seconds. */
    int step = 0;
    int steps_in_minute = ceil_double_to_int(MINUTE / egg->time_step);

    /* Copy initial temperature to saved array. */
    save_data(egg->temp_vector, saved_egg_temps, egg->mesh_size, *saved_data_count);
    *saved_data_count = *saved_data_count + 1;

    while (gsl_vector_get(egg->temp_vector, egg->cooking_point) < TEMP_COOKED) {
        /* Loop until outermost yolk point reaches the temperature at which it is cooked. */
        step += 1;
        egg_setup_rhs_equation(egg);
        error = solve_matrix_equation(egg->matrix, egg->rhs_vector, egg->temp_vector);
        if (error != NO_ERROR) {
            return error;
        }

        if ((step % steps_in_minute) == 0) {
            /* Temperature saved every time a minute passes in simulation. */
            save_data(egg->temp_vector, saved_egg_temps, egg->mesh_size, *saved_data_count);
            *saved_data_count = *saved_data_count + 1;
        }
    }

    /* Copy final temperature to saved array. */
    save_data(egg->temp_vector, saved_egg_temps, egg->mesh_size, *saved_data_count);
    *saved_data_count = *saved_data_count + 1;

    /* Print data of rod temperatures at specific times to a file. */
    error = egg_print_data(egg, saved_egg_temps, *saved_data_count, (step * egg->time_step));
    if (error != NO_ERROR) {
        return error;
    }

    /* Tell user how long it took for the egg to cook. */
    printf("The time taken to cook the egg of mass %.2lfg, was %.3lf seconds.\n", egg->mass, (step * egg->time_step));

    return NO_ERROR;
}

/* Find information on egg and model it. */
Error egg_setup_model(Egg *egg) {
    Error error;

    error = create_egg(egg);
    if (error != NO_ERROR) {
        return error;
    }

    if (egg->hard_boil == true) {
        printf("An egg being hard boiled will be modelled.\n");
    }
    else {
        printf("An egg being soft boiled will be modelled.\n");
    }

    /* Initialise and allocate memory for a large array to store values of temperature later used for plotting. */
    double *saved_egg_temps = malloc(sizeof(double) * egg->mesh_size * EGG_MAX_PLOTS);
    if (saved_egg_temps == NULL) {
        return memory_error_print();
    }

    /* Count to keep a track of number of sets of data saved. */
    int saved_data_count = 0;

    /* Process of simulating the egg being boiled and recording data from it. */
    error = egg_model(egg, saved_egg_temps, &saved_data_count);
    if (error != NO_ERROR) {
        safe_free(saved_egg_temps);
        return error;
    }

    safe_free(saved_egg_temps);

    /* Plot images of simulation for user to see. */
    error = egg_plot(egg, saved_data_count);
    if (error != NO_ERROR) {
        return error;
    }

    printf("Data has been written to the file: %s.\n", egg->file_name);

    return NO_ERROR;
}

/* Set default values for modelling an egg. */
void give_egg_default_values(Egg *egg) {
    egg->hard_boil = false;
    egg->time_step = DEFAULT_EGG_TIME_STEP;
    egg->resolution = DEFAULT_EGG_RESOLUTION;
    egg->mass = DEFAULT_EGG_MASS;
    egg->initial_temperature = ROOM_TEMP;
    strcpy(egg->file_name, DEFAULT_EGG_OUTPUT_FILE);
}


/* Check that each command line value is within acceptable range of values. */
Error check_arg(const char *arg_value, double *struct_value, const double upper_lim, const double lower_lim, const char *arg) {
    Error error;

    error = get_arg_double(arg_value, struct_value, arg);
    if (error != NO_ERROR) {
        return error;
    }
    if (*struct_value > upper_lim || *struct_value <= lower_lim) {
        /* Returns error if not a valid value for rod model. */
        return invalid_args_print("Value given is too large or small.", arg, optarg);
    }

    return NO_ERROR;
}

/* Check that a model has been set before attempting to give other values. */
Error command_line_order_check(const bool *opt_rod, const bool *opt_egg, const char *optarg, const char *arg) {
    /* Check that at least, and only, one option has been set. */
    if ((*opt_rod == false && *opt_egg == false) || (*opt_rod == true && *opt_egg == true)) {
        return invalid_args_print("Please select one model you would like to do and put it first on the command line.", arg, optarg);
    }
    return NO_ERROR;
}

/* Retrieve values from command line arguments, check suitability and give to structures. */
Error get_args(Rod *rod, Egg *egg, const int argc, char *const *argv, bool *opt_rod, bool *opt_egg) {
    Error error;
    int c;

    while (1) {
        static struct option long_options[] = {
                /* These options set a flag. */
                {"verbose", no_argument,       &verbose_flag, 1},
                {"brief",   no_argument,       &verbose_flag, 0},
                /* These options don’t set a flag.
                   We distinguish them by their indices. */
                {"rod",     no_argument,       0, 'n'},
                {"egg",  no_argument,       0, 'e'},
                {"hardboil",  no_argument,       0, 'b'},
                {"timestep",  required_argument, 0, 's'},
                {"resolution",  required_argument, 0, 'r'},
                {"mass",  required_argument, 0, 'm'},
                {"temperature",  required_argument, 0, 't'},
                {"file",  required_argument, 0, 'f'},
                {"help",  no_argument, 0, 'h'},
                {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, ":nebs:r:m:t:f:h", long_options, &option_index);

        /* Detect the end of the options and breaks out of the initial while loop. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf ("\n");
                break;

            case 'n':
                *opt_rod = true;
                give_rod_default_values(rod);
                break;

            case 'e':
                *opt_egg = true;
                give_egg_default_values(egg);
                break;

            case 'b':
                error = command_line_order_check(opt_rod, opt_egg, "", long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                /* Checking for which model is being done. */
                if (*opt_rod == true) {
                    return invalid_args_print("Option selected does not correspond with model chosen.", long_options[option_index].name, "");
                }
                else if (*opt_egg == true) {
                    egg->hard_boil = true;
                }
                break;

            case 's':
                error = command_line_order_check(opt_rod, opt_egg, optarg, long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                /* Checking for which model is being done. */
                if (*opt_rod == true) {
                    error = check_arg(optarg, &rod->time_step, MAX_ROD_TIME_STEP, MIN_TIME_STEP, long_options[option_index].name);
                }
                else if (*opt_egg == true) {
                    error = check_arg(optarg, &egg->time_step, MAX_EGG_TIME_STEP, MIN_TIME_STEP, long_options[option_index].name);
                }
                if (error != NO_ERROR) {
                    return error;
                }
                break;

            case 'r':
                error = command_line_order_check(opt_rod, opt_egg, optarg, long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                /* Checking for which model is being done. */
                if (*opt_rod == true) {
                    error = check_arg(optarg, &rod->resolution, MAX_RESOLUTION, MIN_RESOLUTION, long_options[option_index].name);
                }
                else if (*opt_egg == true) {
                    error = check_arg(optarg, &egg->resolution, MAX_RESOLUTION, MIN_RESOLUTION, long_options[option_index].name);
                }
                if (error != NO_ERROR) {
                    return error;
                }
                break;

            case 'm':
                error = command_line_order_check(opt_rod, opt_egg, optarg, long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                /* Checking for which model is being done. */
                if (*opt_rod == true) {
                    return invalid_args_print("Option selected does not correspond with model chosen.", long_options[option_index].name, optarg);
                }
                else if (*opt_egg == true) {
                    error = check_arg(optarg, &egg->mass, MAX_EGG_MASS, MIN_EGG_MASS, long_options[option_index].name);
                    if (error != NO_ERROR) {
                        return error;
                    }
                }
                break;

            case 't':
                error = command_line_order_check(opt_rod, opt_egg, optarg, long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                /* Checking for which model is being done. */
                if (*opt_rod == true) {
                    return invalid_args_print("Option selected does not correspond with model chosen.", long_options[option_index].name, optarg);
                }
                else if (*opt_egg == true) {
                    error = check_arg(optarg, &egg->initial_temperature, TEMP_COOKED, MIN_EGG_TEMP, long_options[option_index].name);
                    if (error != NO_ERROR) {
                        return error;
                    }
                }
                break;

            case 'f':
                error = command_line_order_check(opt_rod, opt_egg, optarg, long_options[option_index].name);
                if (error != NO_ERROR) {
                    return error;
                }

                if (strlen(optarg) >= FILE_NAME_BUFFER) {
                    /* Checks that the string isn't too long. */
                    return invalid_args_print("File name given is too long.", long_options[option_index].name, optarg);
                }
                else {
                    if (*opt_rod == true) {
                        strcpy(rod->file_name, optarg);
                    }
                    else if (*opt_egg == true) {
                        strcpy(egg->file_name, optarg);
                    }
                }
                break;

            case 'h':
                print_help();
                return HELP_CALLED;

            case ':':
                /* Missing option argument */
                fprintf(stderr, "\nThe option '-%c' requires an argument\n", optopt);
                print_help();
                return INVALID_ARGS;

            case '?':
                /* getopt_long already printed an error message. */
                fprintf(stderr, "Invalid option -%c given.\n", optopt);
                print_help();
                return INVALID_ARGS;

            default:
                return INVALID_ARGS;
        }
    }

    /* Instead of reporting ‘--verbose’ and ‘--brief’ as they are encountered, we report the final status resulting from them. */
    if (verbose_flag)
        puts ("verbose flag is set");
    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        fprintf(stderr, "non-option ARGV-elements: ");
        while (optind < argc)
            fprintf (stderr, "%s ", argv[optind++]);
        fprintf(stderr, "\n");
    }
    return NO_ERROR;
}

/* Free both structures allocated for in main. */
void free_all(Rod *rod, Egg *egg) {
    free_rod(rod);
    free_egg(egg);
}


int main (int argc, char **argv) {
    Error error;
    /* Set counters for use in get_opt to 0. */
    bool opt_rod = false, opt_egg = false;

    /* Allocate memory for structures used in the simulations. */
    Rod *rod = malloc(sizeof(Rod));
    Egg *egg = malloc(sizeof(Egg));
    if (egg == NULL || rod == NULL) {
        free_all(rod, egg);
        return memory_error_print();
    }

    error = get_args(rod, egg, argc, argv, &opt_rod, &opt_egg);
    if (error != NO_ERROR) {
        free_all(rod, egg);
        return error;
    }

    if (opt_rod == true && opt_egg == false) {
        /* If option to model a nuclear waste rod is selected, and no others are, the rod is modelled. */
        error = rod_setup_model(rod);
        if (error != NO_ERROR) {
            free_all(rod, egg);
            return error;
        }
    }
    else if (opt_egg == true && opt_rod == false) {
        /* If option to model an egg is selected, and no others are, the egg is modelled. */
        error = egg_setup_model(egg);
        if (error != NO_ERROR) {
            free_all(rod, egg);
            return error;
        }
    }

    free_all(rod, egg);
    return NO_ERROR;
}

/*
Typical program output:

jemgodden2@Jems-MacBook-Air Thermal_Diffusion % ./compile
jemgodden2@Jems-MacBook-Air Thermal_Diffusion % ./a.out --rod --timestep 0.05 --file olsen_kettle.txt
A nuclear waste rod will be modelled.
Data has been written to the file: olsen_kettle.txt.
jemgodden2@Jems-MacBook-Air Thermal_Diffusion % ./a.out --egg --mass 60.0 --temperature 277.0 --resolution 1.0
An egg being soft boiled will be modelled.
The time taken to cook the egg of mass 60.00g, was 316.300 seconds.
Data has been written to the file: egg_data.txt.
 */