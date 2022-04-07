#!/bin/bash

## Code to run different scripts on any system through an interface rather than direct submission

# Load any necessary modules
module load slurm

#Load the current working directory
export START_DIR=$PWD

# Store variables in case of exit needed from function
trap "exit 1" TERM
export TOP_PID=$$

# Store script location for relative path construction
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Guesses for the root of the repository
export ROOT_DIR_OPT_1="`realpath $SCRIPT_DIR/../`"
export ROOT_DIR_OPT_2=$START_DIR
export ROOT_DIR_OPT_3=$(head -n 1 $SCRIPT_DIR/.cache/stored_dir)

# Select the root directory
echo "Select which of these is the root directory for the  real-time-mcmc repository: (May be multiple, in which case pick any correct directory)"
select ROOT_DIR_SEL in $ROOT_DIR_OPT_1 $ROOT_DIR_OPT_2 $ROOT_DIR_OPT_3 "Other"
do
    case $ROOT_DIR_SEL in

        $ROOT_DIR_OPT_1 | $ROOT_DIR_OPT_2 | $ROOT_DIR_OPT_3)
            echo "Selected $ROOT_DIR_SEL as the root directory of the real-time-mcmc repository."
            export ROOT_DIR=$ROOT_DIR_SEL
            break
            ;;

        "Other")
            echo "Input location of root directory: (Note use full path i.e. ~ expansion not yet working)"
            read ROOT_DIR
            case $ROOT_DIR in

                ?*)
                    if [ ! -d $ROOT_DIR ] 
                    then
                        echo "Directory not found. Exiting..."
                        exit 0
                    fi
                    ROOT_DIR="`realpath $ROOT_DIR`"
                    echo "Using the following option: $ROOT_DIR"
                    options="$options --array=$array"
                    echo "Would you like to replace option 3 ($ROOT_DIR_OPT_3) with $ROOT_DIR in future?"
                    select save_root_opt in Yes No
                    do
                        case $save_root_opt in

                            Yes)
                                echo "Saving $ROOT_DIR..."
                                echo $ROOT_DIR >| $SCRIPT_DIR/.cache/stored_dir
                                echo "Saved"
                                break
                                ;;

                            No)
                                echo "Chose not to save option. Continuing..."
                                break
                                ;;

                        esac
                    done
                    break
                    ;;

                *)
                    echo "No directory input. Exiting..."
                    exit 0
                    ;;

            esac
            break
            ;;


    esac
done

# Pick a script to run
echo "Select an option for a script to run (or select 6 to exit):"

select script_opt in "Pre-process" "Run Model" "Post-process" "Save Endstates" "Restart Runs" "Exit to Shell"
do
    case $script_opt in

        "Pre-process")
            export RUN_NAME="Pre-process inputs"
            export SLURM_LOC=$ROOT_DIR/submission_scripts/submit_cu_hpc_input_process
            break
            ;;
        
        "Run Model")
            export RUN_NAME="Run model"
            export SLURM_LOC=$ROOT_DIR/submission_scripts/submit_cu_hpc
            break
            ;;
        
        "Post-process")
            export RUN_NAME="Post-process outputs"
            export SLURM_LOC=$ROOT_DIR/submission_scripts/submit_post_process
            break
            ;;
        
        "Save Endstates")
            export RUN_NAME="Save Endstates"
            export SLURM_LOC=$ROOT_DIR/submission_scripts/submit_save_endstates
            break
            ;;
        
        "Restart Runs")
            export RUN_NAME="Restart Runs"
            export SLURM_LOC=$ROOT_DIR/submission_scripts/submit_restart
            break
            ;;
        
        "Exit to Shell")
            echo "Exiting to Shell"
            cd $START_DIR
            exit 0
            ;;

        *) 
            echo "Unrecognised selection. Try again with another option"
            ;;
    esac
done

export SLURM_LOC="`realpath $SLURM_LOC`"


# Confirm that the script chosen is the one wanted
echo "Please confirm that you want to run the $RUN_NAME script at $SLURM_LOC:"

select opt_confirm in Yes No
do
    case $opt_confirm in
        
        Y | y | Yes | yes)
            echo "Confirmed selection of $RUN_NAME and continuing to selection of input parameters."
            break
            ;;

        N | n | No | no)
            echo "Cancelled selection. Re-run script with correct selection."
            cd $START_DIR
            exit 0
            ;;
        *)
            echo "Unrecognised command. Try again with another option"
            ;;
    esac
done

export options=""

# For array based scripts set the array variable
case $script_opt in

    "Run Model" | "Post-process" | "Save Endstates" | "Restart Runs")
        echo "What should be used as the array inputs? (press Enter to use default args from file)"
        read array
        case $array in

            ?*)
                echo "Using the following option: --array=$array"
                options="$options --array=$array"
                ;;

            *?)
                "Using Default args from file"
                ;;

        esac
        ;;

esac

# Set working directory
echo "Where should the script be run from (i.e. for run model needs to be the folder {root}/model_runs/{date})?
Note: If nothing sent uses current directory $START_DIR
Start from root by starting directory with /
Otherwise starts from relative path $ROOT_DIR"
read CWDIR
case $CWDIR in

    ?*)
        firstCharacter=${CWDIR:0:1}
        if [[ ! firstCharacter == "/" ]]
        then
            CWDIR="`realpath $ROOT_DIR/$CWDIR`"
        fi
        echo $CWDIR
        if [ ! -d $CWDIR ] 
        then
            echo "Directory not found. Exiting..."
            exit 0
        fi
        CWDIR="`realpath $CWDIR`"
        echo "Using the following as the directory: $CWDIR"
        cd $CWDIR
        ;;

    *)
        echo "Running from $START_DIR"
        cd $START_DIR
        ;;

esac

# Run the commands
export CMD="sbatch $options $SLURM_LOC"

echo "CMD = $CMD"

eval $CMD

# Return to starting directory
cd $START_DIR