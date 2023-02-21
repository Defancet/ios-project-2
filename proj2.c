#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/mman.h>
#include <sys/types.h>

typedef struct args {    // struct to hold arguments
    int NO;
    int NH;
    int TI;
    int TB;
} args_t;

typedef struct semaphores {    // struct to hold semaphores
    sem_t hyd;
    sem_t oxy;
    sem_t mol;
    sem_t mutex;
    sem_t mutex2;
    sem_t br;
    sem_t br2;
    sem_t out;
    sem_t end;
} semaphores_t;

typedef struct shared_memory {    // struct to hold shared memory
    int num_hydrogen;
    int num_oxygen;
    int curr_hydrogen;
    int curr_oxygen;
    int atom_id;
    int mol_id;
    int cnt;
    int next;
} shared_memory_t;

FILE *file;    // file pointer

void print_to_file(int state, semaphores_t *semaphores, shared_memory_t *shared_memory, char atom, int id);    // function to print to file
int args_process(int argc, char **argv, args_t *args);    // function to process arguments
void create_atom(semaphores_t *semaphores, shared_memory_t *shared_memory);    // function to create atom
void next(semaphores_t *semaphores, shared_memory_t *shared_memory);    // function to move atom to next position
void stop(semaphores_t *semaphores, shared_memory_t *shared_memory);    // function to prevent atom from passing through
int hydrogen(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args, int id);    // function to create hydrogen
int oxygen(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args, int id);    // function to create oxygen
void destroy_sem(semaphores_t *semaphores);    // function to destroy semaphores
int init(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args);    // function to initialize shared memory and semaphores

int main(int argc, char **argv) {
    args_t args;    // struct to hold arguments
    pid_t pid;    // pid of child process
    shared_memory_t *shared_memory = mmap(NULL, sizeof(shared_memory_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
    semaphores_t *semaphores = mmap(NULL, sizeof(shared_memory_t), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, 0, 0);
    if (shared_memory == MAP_FAILED) {    // check if mmap failed to allocate shared memory
        fprintf(stderr, "Unable to create shared memory.");
        return 1;
    }

    if (args_process(argc, argv, &args) == 1) {    // check if arguments were processed correctly
        return 1;
    }

    if ((file = fopen("proj2.out", "w")) == NULL) {    // open file to write to and check if it was opened correctly
        fprintf(stderr, "Error when opening a file");
        return 1;
    }

    init(semaphores, shared_memory, args); // initialize shared memory and semaphores

    for (int i = 0; i < args.NH; i++) {    // create hydrogen atoms
        pid = fork();    // create child process
        if (pid == -1) {    // check if fork failed to create child process
            fprintf(stderr, "Fork error\n");
            return 1;
        }
        if (pid == 0) {    // child process
            hydrogen(semaphores, shared_memory, args, i + 1);    // create hydrogen atom
            return 1;
        }
    }

    for (int i = 0; i < args.NO; i++) {    // create oxygen atoms
        pid = fork();    // create child process
        if (pid == -1) {    // check if fork failed to create child process
            fprintf(stderr, "Fork error\n");
            return 1;
        }
        if (pid == 0) {     // child process
            oxygen(semaphores, shared_memory, args, i + 1);     // create oxygen atom
            return 1;
        }
    }

    while (wait(NULL) > 0);    // wait for all child processes to finish

    destroy_sem(semaphores);    // destroy semaphores
    fclose(file);    // close file

    if (munmap(shared_memory, sizeof(shared_memory_t)) == -1) {    // unmap shared memory and check if it was unmapped correctly
        fprintf(stderr, "Unable to unmap shared memory.");
        return 1;
    }
    if (munmap(semaphores, sizeof(semaphores_t)) == -1) {    // unmap semaphores and check if it was unmapped correctly
        fprintf(stderr, "Unable to unmap semaphores.");
        return 1;
    }
    return 0;
}

void print_to_file(int state, semaphores_t *semaphores, shared_memory_t *shared_memory, char atom, int id) {    // function to print to file
    sem_wait(&semaphores->out);    // wait for access to output file
        if (state == 1) {
            fprintf(file, "%d: %c %d: started\n", ++shared_memory->atom_id, atom, id);    // print to file that atom is started
        }
        else if (state == 2) {
            fprintf(file, "%d: %c %d: going to queue\n", ++shared_memory->atom_id, atom, id);    // print to file that atom is going to queue
        }
        else if (state == 3) {
            fprintf(file, "%d: %c %d: creating molecule %d\n", ++shared_memory->atom_id, atom, id, shared_memory->mol_id / 3);    // print to file that molecule is being created
        }
        else if (state == 4) {
            fprintf(file, "%d: %c %d: molecule %d created\n", ++shared_memory->atom_id, atom, id, shared_memory->mol_id++ / 3);    // print to file that molecule is created
        }
        else if (state == 5) {
            fprintf(file, "%d: %c %d: not enough O or H\n", ++shared_memory->atom_id, atom, id);    // print to file that there is not enough atoms of O or H
        }
        else if (state == 6) {
            fprintf(file, "%d: %c %d: not enough H\n", ++shared_memory->atom_id, atom, id);    // print to file that there is not enough atoms of H
        }
    fflush(file);    // flush file
    sem_post(&semaphores->out);    // release access to output file
}

int args_process(int argc, char **argv, args_t *args) {    // function to process arguments
    char *tmp_no = "", *tmp_nh = "", *ptr_ti = "", *ptr_tb = "";    // temporary variables
    if (argc != 5) {    // check if there are 4 arguments
        fprintf(stderr, "Invalid number of arguments\n");
        return 1;
    }
    args->NO = (int)strtol(argv[1], &tmp_no, 10);    // number of oxygen molecules
    args->NH = (int)strtol(argv[2], &tmp_nh, 10);    // number of hydrogen molecules
    args->TI = (int)strtol(argv[3], &ptr_ti, 10);    // time interval for hydrogen molecules
    args->TB = (int)strtol(argv[4], &ptr_tb, 10);    // time interval for oxygen molecules

    if ((strcmp(tmp_no, "") != 0) || args->NO <= 0) {    // check if number of oxygen molecules is valid
        fprintf(stderr, "Incorrect NO argument\n");
        return 1;
    }
    if ((strcmp(tmp_nh, "") != 0) || args->NH <= 0) {    // check if number of hydrogen molecules is valid
        fprintf(stderr, "Incorrect NH argument\n");
        return 1;
    }
    if ((strcmp(ptr_ti, "") != 0) || args->TI < 0 || args->TI > 1000) {    // check if time interval for hydrogen molecules is valid
        fprintf(stderr, "Incorrect TI argument\n");
        return 1;
    }
    if ((strcmp(ptr_tb, "") != 0) || args->TB < 0 || args->TB > 1000) {    // check if time interval for oxygen molecules is valid
        fprintf(stderr, "Incorrect TB argument\n");
        return 1;
    }
    return 0;
}

void create_atom(semaphores_t *semaphores, shared_memory_t *shared_memory) {    // function to create atom
    sem_post(&(semaphores->hyd));
    sem_post(&(semaphores->hyd));
    sem_post(&(semaphores->oxy));
    shared_memory->num_hydrogen -= 2;
    shared_memory->curr_hydrogen -= 2;
    shared_memory->num_oxygen -= 1;
    shared_memory->curr_oxygen -= 1;
}

void next(semaphores_t *semaphores, shared_memory_t *shared_memory) {    // function to move atom to next position
    if (shared_memory->curr_hydrogen < 2 || shared_memory->curr_oxygen < 1) {    // check if there are enough atoms to move
        shared_memory->next = 1;    // set next to 1 to indicate that there are no more atoms to move
        for (int i = 0; i < shared_memory->curr_oxygen; i++) {    // loop to move all oxygen atoms
            sem_post(&(semaphores->oxy));
            sem_post(&(semaphores->mutex));
        }
        for (int i = 0; i < shared_memory->curr_hydrogen; i++) {    // loop to move all hydrogen atoms
            sem_post(&(semaphores->hyd));
            sem_post(&(semaphores->mutex));
        }
    }
}

void stop(semaphores_t *semaphores, shared_memory_t *shared_memory) {    // function to prevent atoms from passing through
    sem_wait(&(semaphores->mutex2));
    shared_memory->cnt += 1;
    if (shared_memory->cnt == 3) {    // check if all atoms are created
        for (int i = 0; i < 3; i++)
            sem_post(&(semaphores->br));
    }
    sem_post(&(semaphores->mutex2));
    sem_wait(&(semaphores->br));
    sem_wait(&(semaphores->mutex2));
    shared_memory->cnt -= 1;
    if (shared_memory->cnt == 0) {    // check if all atoms are created
        for (int i = 0; i < 3; i++)
            sem_post(&(semaphores->br2));
    }
    sem_post(&(semaphores->mutex2));
    sem_wait(&(semaphores->br2));
}

int hydrogen(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args, int id) {    // function to create hydrogen
    print_to_file(1, semaphores, shared_memory, 'H', id);    // call the function to print to file that atom of hydrogen is started

    usleep((rand() % (args.TI + 1)) * 1000);    // sleep for random time interval

    print_to_file(2, semaphores, shared_memory, 'H', id);    // call the function to print to file that atom of hydrogen is going to queue

    sem_wait(&(semaphores->mutex));
    if (shared_memory->next) {    // check if there are no more atoms to move
        print_to_file(5, semaphores, shared_memory, 'H', id);    // call the function to print to file that there is not enough atoms of O or H

        shared_memory->curr_hydrogen--;    // decrease number of hydrogen atoms

        sem_post(&(semaphores->mutex));
        sem_post(&(semaphores->hyd));
        sem_post(&(semaphores->oxy));
        return 1;
    }

    if (shared_memory->curr_hydrogen < 2 || shared_memory->curr_oxygen < 1) {    // check if there are enough hydrogen and oxygen atoms to create molecule
        print_to_file(5, semaphores, shared_memory, 'H', id);
        shared_memory->curr_hydrogen--;
        sem_post(&(semaphores->mutex));
        return 1;
    }

    shared_memory->num_hydrogen++;

    if (shared_memory->num_hydrogen >= 2 && shared_memory->num_oxygen >= 1) {    // check if there are enough hydrogen and oxygen atoms to create atom
        create_atom(semaphores, shared_memory);
    } else {    // if there are not enough hydrogen and oxygen atoms to create atom
        sem_post(&(semaphores->mutex));
    }

    sem_wait(&(semaphores->hyd));
    if (shared_memory->next) {    // check if there are no more atoms to move
        print_to_file(5, semaphores, shared_memory, 'H', id);    // call the function to print to file that there is not enough atoms of O or H

        shared_memory->curr_hydrogen--;    // decrease number of hydrogen atoms

        sem_post(&(semaphores->mutex));
        sem_post(&(semaphores->hyd));
        sem_post(&(semaphores->oxy));
        return 1;
    }

    print_to_file(3, semaphores, shared_memory, 'H', id);    // call the function to print to file that molecule is being created

    stop(semaphores, shared_memory);    // call the function to prevent atom from going through

    print_to_file(4, semaphores, shared_memory, 'H', id);    // call the function to print to file that molecule H is created

    sem_post(&(semaphores->end));
    return 0;
}

int oxygen(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args, int id) {    // function to create oxygen
    print_to_file(1, semaphores, shared_memory, 'O', id);    // call the function to print to file that oxygen atom is started

    usleep((rand() % (args.TI + 1)) * 1000);    // sleep for random time interval

    print_to_file(2, semaphores, shared_memory, 'O', id);    // call the function to print to file that atom of oxygen is going to queue

    sem_wait(&(semaphores->mutex));
    if (shared_memory->next) {    // check if there are no more atoms to move
        print_to_file(6, semaphores, shared_memory, 'O', id);    // call the function to print to file that there is not enough atoms of H

        shared_memory->curr_oxygen--;    // decrease number of oxygen atoms

        sem_post(&(semaphores->mutex));
        sem_post(&(semaphores->hyd));
        sem_post(&(semaphores->oxy));
        return 1;
    }

    if (shared_memory->curr_hydrogen < 2 || shared_memory->curr_oxygen < 1) {    // check if there are enough hydrogen and oxygen atoms to create molecule
        print_to_file(6, semaphores, shared_memory, 'O', id);
        shared_memory->curr_oxygen--;
        sem_post(&(semaphores->mutex));
        return 1;
    }

    shared_memory->num_oxygen++;    // increase number of oxygen atoms
    if (shared_memory->num_hydrogen >= 2 && shared_memory->num_oxygen >= 1) {    // check if there are enough hydrogen atoms to create atom
        create_atom(semaphores, shared_memory);
    } else {    // if there are not enough hydrogen atoms to create atom
        sem_post(&(semaphores->mutex));
    }

    sem_wait(&(semaphores->oxy));
    if (shared_memory->next) {    // check if there are no more atoms to move
        print_to_file(6, semaphores, shared_memory, 'O', id);    // call the function to print to file that there is not enough atoms of H

        shared_memory->curr_oxygen--;    // decrease number of oxygen atoms

        sem_post(&(semaphores->mutex));
        sem_post(&(semaphores->hyd));
        sem_post(&(semaphores->oxy));
        return 1;
    }

    print_to_file(3, semaphores, shared_memory, 'O', id);    // call the function to print to file that molecule O is being created

    usleep((rand() % (args.TB + 1)) * 1000);    // sleep for random time interval

    stop(semaphores, shared_memory);    // call the function to prevent atom from going through

    sem_wait(&(semaphores->end));
    sem_wait(&(semaphores->end));

    print_to_file(4, semaphores, shared_memory, 'O', id);    // call the function to print to file that molecule O is created

    next(semaphores, shared_memory);    // call the function to check if there are more atoms to move

    sem_post(&(semaphores->mutex));
    return 0;
}

void destroy_sem(semaphores_t *semaphores) {    // function to destroy semaphores
    sem_destroy(&(semaphores->hyd));
    sem_destroy(&(semaphores->oxy));
    sem_destroy(&(semaphores->mol));
    sem_destroy(&(semaphores->mutex));
    sem_destroy(&(semaphores->mutex2));
    sem_destroy(&(semaphores->br));
    sem_destroy(&(semaphores->br2));
    sem_destroy(&(semaphores->out));
    sem_destroy(&(semaphores->end));
}

int init(semaphores_t *semaphores, shared_memory_t *shared_memory, args_t args) {    // function to initialize shared memory and semaphores
    sem_init(&(semaphores->oxy), 1, 0);    // initialize semaphores
    sem_init(&(semaphores->hyd), 1, 0);
    sem_init(&(semaphores->mol), 1, 0);
    sem_init(&(semaphores->mutex), 1, 1);
    sem_init(&(semaphores->mutex2), 1, 1);
    sem_init(&(semaphores->br), 1, 0);
    sem_init(&(semaphores->br2), 1, 0);
    sem_init(&(semaphores->out), 1, 1);
    sem_init(&(semaphores->end), 1, 0);

    shared_memory->num_hydrogen = 0;    // initialize shared memory
    shared_memory->num_oxygen = 0;
    shared_memory->curr_hydrogen = args.NH;
    shared_memory->curr_oxygen = args.NO;
    shared_memory->atom_id = 0;
    shared_memory->mol_id = 3;
    shared_memory->cnt = 0;
    shared_memory->next = 0;

    return 0;
}
