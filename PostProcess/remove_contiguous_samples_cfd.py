#!/usr/bin/env python
"""
Created on Fri Aug 28 15:58:04 2020

@author: francesco
"""
import sys
import time
import shutil


def get_best_targets(
    cluster,
    best_list,
    discard_list,
    tau,
    chrom,
    pos,
    total,
    true_guide,
    snp_info,
    cfd,
    sort_order,
):
    # avoid crush when cluster is empty in the first call
    if not cluster:
        return

    list_ref = list()
    dict_var = dict()
    for ele in cluster:
        if ele[snp_info] == "n":
            list_ref.append(ele)
        else:
            # merge samples of identical targets (coming from different VCF datasets)
            if (ele[pos], ele[snp_info]) in dict_var.keys():
                dict_var[(ele[pos], ele[snp_info])][0][true_guide - 2] = (
                    dict_var[(ele[pos], ele[snp_info])][0][true_guide - 2]
                    + ","
                    + ele[true_guide - 2]
                )  # true_guide - 2 points to samples column
                dict_var[(ele[pos], ele[snp_info])][0][snp_info - 2] = (
                    dict_var[(ele[pos], ele[snp_info])][0][snp_info - 2]
                    + ","
                    + ele[snp_info - 2]
                )  # snp_info - 2 points to rsID column
                dict_var[(ele[pos], ele[snp_info])][0][snp_info - 1] = (
                    dict_var[(ele[pos], ele[snp_info])][0][snp_info - 1]
                    + ","
                    + ele[snp_info - 1]
                )  # ttuesnp_info_guide - 2 points to AF column
            else:
                dict_var[(ele[pos], ele[snp_info])] = [ele]

    final_list_best_ref = list()
    var_only = False
    for target in list_ref:
        final_list_best_ref.append(target)
    if not final_list_best_ref:
        var_only = True

    final_list_best_var = list()
    # for each snp_info in dict, extract the targets
    for key in dict_var.keys():
        list_var = dict_var[key]
        # copy the targets in the variant list, adding unique if no ref target is found
        for target in list_var:
            if var_only:
                target[12] = "y"
            final_list_best_var.append(target)

    temp_final_list_best_var = list()
    # for target in final_list_best_var:
    for target in final_list_best_var:
        # remove duplicates into snp info col
        target[snp_info] = ",".join(set(target[snp_info].split(",")))
        # remove duplicate into rsID col
        target[snp_info - 2] = ",".join(set(target[snp_info - 2].split(",")))
        # remove duplicate into AF col
        target[snp_info - 1] = ",".join(set(target[snp_info - 1].split(",")))
        # remove duplicate into samples col
        target[true_guide - 2] = ",".join(set(target[true_guide - 2].split(",")))
        # append to temp list
        temp_final_list_best_var.append(target)

    # final list with polished targets (no duplicates in snp data)
    final_list_best_var = temp_final_list_best_var

    # check if lists are empty
    validity_check_ref = False
    validity_check_var = False
    if final_list_best_ref:
        validity_check_ref = True
    if final_list_best_var:
        validity_check_var = True

    # extract best target for each criteria
    if sort_order == "score":
        # sort per score (CFD or CRISTA)
        if validity_check_ref:
            final_list_best_ref = sorted(
                final_list_best_ref, key=lambda x: (-float(x[cfd]), int(x[total]))
            )
        if validity_check_var:
            final_list_best_var = sorted(
                final_list_best_var, key=lambda x: (-float(x[cfd]), int(x[total]))
            )
        if var_only:  # no ref found
            # count the residual targets in the list
            final_list_best_var[0][cfd - 1] = str(len(final_list_best_var) - 1)
            # append the best target to best_file
            best_list.append(final_list_best_var[0])
            # pop the best target from the list
            bestTarget = final_list_best_var.pop(0)
        elif validity_check_ref and validity_check_var:  # ref and var targets found
            if float(final_list_best_ref[0][cfd]) >= float(final_list_best_var[0][cfd]):
                final_list_best_ref[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                best_list.append(final_list_best_ref[0])
                bestTarget = final_list_best_ref.pop(0)
            else:
                final_list_best_var[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                best_list.append(final_list_best_var[0])
                bestTarget = final_list_best_var.pop(0)
        else:  # only ref
            final_list_best_ref[0][cfd - 1] = str(len(final_list_best_ref) - 1)
            best_list.append(final_list_best_ref[0])
            bestTarget = final_list_best_ref.pop(0)
        # write all the remaining targets in the alt file
        for count, elem in enumerate(final_list_best_ref):
            final_list_best_ref[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            discard_list.append(elem)
        for count, elem in enumerate(final_list_best_var):
            final_list_best_var[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            discard_list.append(elem)
    else:
        # sort for total (mm+bul) in target
        if validity_check_ref:
            final_list_best_ref = sorted(
                final_list_best_ref,
                key=lambda x: (int(x[total - 2]), int(x[total - 1])),
            )
        if validity_check_var:
            final_list_best_var = sorted(
                final_list_best_var,
                key=lambda x: (int(x[total - 2]), int(x[total - 1])),
            )
        if var_only:  # no ref found
            # count the residual targets in the list
            final_list_best_var[0][cfd - 1] = str(len(final_list_best_var) - 1)
            # append the best target to best_file
            best_list.write(final_list_best_var[0])
            # pop the best target from the list
            bestTarget = final_list_best_var.pop(0)
        elif validity_check_ref and validity_check_var:  # ref and var targets found
            if int(final_list_best_ref[0][total]) <= int(final_list_best_var[0][total]):
                final_list_best_ref[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                best_list.append(final_list_best_ref[0])
                bestTarget = final_list_best_ref.pop(0)
            else:
                final_list_best_var[0][cfd - 1] = str(
                    len(final_list_best_ref) + len(final_list_best_var) - 1
                )
                best_list.append(final_list_best_var[0])
                bestTarget = final_list_best_var.pop(0)
        else:  # only ref
            final_list_best_ref[0][cfd - 1] = str(len(final_list_best_ref) - 1)
            best_list.append(final_list_best_ref[0])
            bestTarget = final_list_best_ref.pop(0)
        # write all the remaining targets in the alt file
        for count, elem in enumerate(final_list_best_ref):
            final_list_best_ref[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            best_list.append(elem)
        for count, elem in enumerate(final_list_best_var):
            final_list_best_var[count][cfd - 1] = str(
                len(final_list_best_ref) + len(final_list_best_var) - 1
            )
            best_list.append(elem)

    return best_list, discard_list


def merge_results(target_list: list, tau: int, sort_order: str) -> list:
    best_list_final = list()
    discard_list_final = list()
    tmp_best_list = list()
    tmp_discard_list = list()

    chrom = 4  # chromosome position in target
    pos = 6  # position of target
    total = 10  # mm+bul value position
    true_guide = 15  # real guide used in the search
    snp_info = 18  # snp_info position(ref_alt_allele)
    cfd = 20  # score position

    prev_pos = -(tau + 1)
    best_row = ""
    prev_guide = ""
    prev_chr = ""
    prev_snp = ""
    cluster = list()

    for line in target_list:
        splitted = line
        print(splitted)
        if (
            prev_guide != splitted[true_guide]
            or prev_chr != splitted[chrom]
            or int(splitted[pos]) - prev_pos > tau
        ):
            tmp_best_list, tmp_discard_list = get_best_targets(
                cluster,
                list(),
                list(),
                tau,
                chrom,
                pos,
                total,
                true_guide,
                snp_info,
                cfd,
                sort_order,
            )  # type: ignore

            cluster = [splitted]
        else:
            cluster.append(splitted)
        prev_guide = splitted[true_guide]
        prev_pos = int(splitted[pos])
        prev_chr = splitted[chrom]
        prev_snp = splitted[snp_info]
        best_list_final.extend(tmp_best_list)
        discard_list_final.extend(tmp_discard_list)

    tmp_best_list, tmp_discard_list = get_best_targets(
        cluster,
        list(),
        list(),
        tau,
        chrom,
        pos,
        total,
        true_guide,
        snp_info,
        cfd,
        sort_order,
    )  # type: ignore

    best_list_final.extend(tmp_best_list)
    discard_list_final.extend(tmp_discard_list)

    return [best_list_final, discard_list_final]


# start = time.time()
# with open(sys.argv[1], "r") as fileIn:
#     header = fileIn.readline()
#     with open(sys.argv[2], "w") as fileOut:
#         with open(sys.argv[2] + ".discarded_samples", "w") as fileOut_disc:
#             fileOut.write(header)
#             fileOut_disc.write(header)


# shutil.move(sys.argv[2], sys.argv[1])## uncomment to overwrite the input file
# print("Merging done in: " + str(time.time() - start))
