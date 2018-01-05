#!/bin/sh

out_is_correct () { #out_is_correct src_prefix out_prefix test_id [interpop]
    POP_LIST=$(cut -f 3 ${1}.sample | sort | uniq)
    for i in ${POP_LIST}
    do
	[ -e ${2}_${i}.out ] && {
	    cmp -s ${2}_${i}.out ${3}_${i}.ref || {
		echo "File ${3}_${i}.ref has changed" && exit 1
	    }
	}
    done
    [ -n "${4}" ] && {
        { cmp -s ${2}.interpop ${3}.interpop.ref && rm ${2}.interpop ; } || {
	    echo "File ${3}.interpop has changed" && exit 1
	}
    }
    for i in ${POP_LIST} ; do rm ${2}_${i}.out ${2}_${i}_excluded.out ; done
}

make_ref () { #make_ref src_prefix out_prefix test_id [interpop]
    POP_LIST=$(cut -f 3 ${1}.sample | sort | uniq)
    for i in ${POP_LIST}
    do
	mv ${2}_${i}.out ${3}_${i}.ref
    done
    [ -n "${4}" ] && mv ${2}.interpop ${3}.interpop.ref || exit 0
}


regenerate () { #regenerate ops src_prefix out_prefix test_id [interpop]
    ../src/selink ${1} -o "${3}" "${2}" 2>/dev/null || exit 1
    make_ref "${2}" "${3}" ${4} ${5}
}

check () { #check ops src_prefix out_prefix test_id [interpop]
    ../src/selink ${1} -o "${3}" "${2}" 2>/dev/null || exit 1
    out_is_correct "${2}" "${3}" ${4} ${5}
}
