#!/usr/bin/env bash
###############################################################################
#
# Scritp to correctly launch MELT in a docker container.
#
###############################################################################



# Arguments validation.
if (( $# < 4 )); then
	echo -e "Error in $0: Not enough arguments." >&2
	exit 1;
fi


# Configuration file.
CONFIG_FILE=config.yml

# Reference paths.
IMG_BASE=/home

IMG_CONFIG=$IMG_BASE/config
IMG_INPUTS=$IMG_BASE/inputs
IMG_OUTPUTS=$IMG_BASE/outputs
IMG_ASSETS=$IMG_BASE/assets
IMG_REFERENCE=$IMG_BASE/reference
IMG_ANNOTATION=$IMG_BASE/annotation
IMG_TMP=$IMG_BASE/tmp
IMG_EXTRA=



###############################################################################

# Declare arrays to store flags' values.
declare -A PATH_H DIR_H FILE_H

# Get essential values.
while getopts ":c:i:o:d:r:a:t:" OPTION; do
	case "$OPTION" in
		c) PATH_H["c"]=$OPTARG ;;
		i) PATH_H["i"]=$OPTARG ;;
		o) PATH_H["o"]=$OPTARG ;;
		d) PATH_H["d"]=$OPTARG ;;
		r) PATH_H["r"]=$OPTARG ;;
		a) PATH_H["a"]=$OPTARG ;;
		t) PATH_H["t"]=$OPTARG ;;
	esac
done

# Verify values.
for var in $(echo -e "c\ti\to\td\tr\ta\tt"); do
	if [[ ! -v PATH_H[$var] ]]; then
		echo -e "No flag '$var' defined." >&2
		exit 1
	fi
done

# Transform values.
for var in $(echo -e "c\ti\to\td\tr\ta\tt"); do
	path="${PATH_H[$var]}"
	DIR_H[$var]="${path%/*}"
	FILE_H[$var]="${path##*/}"
	
	# Debug code.
	#echo "dir ${DIR_H[$var]} file ${FILE_H[$var]}"
done

# Write config file for make consumption.
cat << EOF > $CONFIG_FILE
EOF



# Run docker!
docker run \
	--rm -ti \
	-u $(id -u):$(id -g) \
	-v "${DIR_H["c"]}:$IMG_CONFIG" \
	-v "${DIR_H["i"]}:$IMG_INPUTS" \
	-v "${DIR_H["o"]}:$IMG_OUTPUTS" \
	-v "${DIR_H["d"]}:$IMG_ASSETS" \
	-v "${DIR_H["r"]}:$IMG_REFERENCE" \
	-v "${DIR_H["a"]}:$IMG_ANNOTATION" \
	-v "${DIR_H["t"]}:$IMG_TMP" \
	melt
#EOF
