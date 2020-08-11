import gstree as gs


def all_valid_delns(seq, leftnick, rightnick, max_size, min_size=1):
    """
    Identifies the complete set of deletions in the
    specified size range that could be formed by:
    *deleting bases from the DNA ends and then ss or ds
    blunt ligating.  this corresponds to all deletions
    spanning or contiguous with a breakpoint.
    *MMEJ between overhangs resulting in a deletion that
    is between the overhangs but not contiguous to either.
    Assumes that 2 overhangs are on opposite strands and
    that overhangs CANNOT be ligated without MMEJ.
    Takes as arguments the top strand of a DNA sequence
    (5'-3') and two slice indices representing positions
    of nicks.  If nicks are the same, the DSB is blunt ended.
    If not, the DSB has overhangs.
    Does not distinguish 3' and 5' overhangs.
    Returns a list of tuples of slice indices
    indicating deletion boundaries.
    seq, leftnick, rightnick, maxSize ="CCCTGTTATCCCTA", 5, 9, 1
    print AllValidDelns(seq, leftnick, rightnick, maxSize, minSize=1)
    [(4, 5), (5, 6), (6, 7), (8, 9), (9, 10)]
        """
    assert leftnick<=rightnick #leftnick must be smaller than or equal to rightnick
    allDelns=[]
    for deln in range(min_size, max_size+1):
        leftWinEdge=leftnick-deln
        rightWinEdge=leftnick
        #window now demarcates leftmost possible deletion of size deln
        while leftWinEdge<=rightnick:
            #if it is not true that the window is between the 2 nicks and touching neither
            if not(leftWinEdge>leftnick and rightWinEdge<rightnick):
                #then add the deletion to the list
                allDelns.append((leftWinEdge, rightWinEdge))
            #if the deletion is between the breakpoints and touching
            #neither, it is only valid if it is a microhomology,
            #ie if the sequence in the deletion window is repeated
            #anywhere else in the overhangs (L,R)
            if (leftWinEdge>leftnick and rightWinEdge<rightnick):
                L=seq[leftnick:leftWinEdge]
                R=seq[rightWinEdge:rightnick]
                if seq[leftWinEdge:rightWinEdge] in L or seq[leftWinEdge:rightWinEdge] in R:
                    allDelns.append((leftWinEdge, rightWinEdge))
             #slide the window 1 base to the right
            leftWinEdge+=1 
            rightWinEdge+=1
    return allDelns


def unique_maxpairs(seq, ML=2):
    """
    Interface with GSTree
    Takes a DNA sequence (string) and returns a list of maximal pairs
    denoted by a tuple of form (index1, index2, wordlength, orientation).
    index1 and index2 are the leftmost indices of the repeated motifs.
    they will be sorted such that index1<=index2.
    wordlength is the length of the repeated motifs.
    orientation is a string, either "fwd" or "rc."
    Optional argument "ML" is minimum length of repeated motif.
    Duplicate pairs (produced because the suffix tree finds rev compls
    by looking for direct repeats between seq and its rev compl) are
    filtered out, but the degenerate case of all sequences being their
    own repeats is not.  Pairs that are both forward and rev compl will
    be returned twice, once with each orientation tag.
    """
    #make the suffix tree t with seq and its rev compl
    t = gs.GSTree()
    t.addstring("fwd", seq)
    rcseq = rev_compl(seq)
    t.addstring("rc", rcseq)
    allmaxpairs = t.maxpairs(minlen=ML)
    #get rid of the rc/fwd and rc/rc maxpairs (redundant)
    #convert maxpairs to form (index1, index2, wordlen, orientation) for easier sorting later
    filtered_maxpairs = []
    for maxpair in allmaxpairs:
        key1, key2 = maxpair[0][0], maxpair[1][0]
        if key1 == "fwd" and key2 == "fwd":
            new = (maxpair[0][1], maxpair[1][1], maxpair[2], "fwd")
            filtered_maxpairs.append(new)
        if key1 == "fwd" and key2 == "rc":
            #first convert both indices to be relative to the forward string
            converted_maxpair = rc2fwd(maxpair, len(seq))
            new = (converted_maxpair[0][1], converted_maxpair[1][1], converted_maxpair[2], "rc")
            filtered_maxpairs.append(new)
    #internally sort the indices to be in ascending order
    sorted_maxpairs = []
    for maxpair in filtered_maxpairs:
        if maxpair[0] > maxpair[1]:
            sorted_maxpair = (maxpair[1], maxpair[0], maxpair[2], maxpair[3])
            sorted_maxpairs.append(sorted_maxpair)
        else: sorted_maxpairs.append(maxpair)
    #remove the duplicate maxpairs
    unique_maxpairs = list(set(sorted_maxpairs))
    #sort the maxpairs relative to each other by asceding left motif edge
    sorted_unique = sorted(unique_maxpairs, key=lambda maxpair:maxpair[0])
    return sorted_unique


def rc2fwd(maxpair, seqlen):
    """
    helper function for unique_maxpairs
    maxpair is as output by GSTree
    converts the index of the second motif in the
    maxpair to be relative to the rc of what it was originally
    relative to and returns the converted maxpair
    """
    rc_index = maxpair[1][1]
    wordlen = maxpair[2]
    fwd_index = seqlen - (rc_index + wordlen)
    new_maxpair = ((maxpair[0][0], maxpair[0][1]), ("fwd", fwd_index), maxpair[2])
    return new_maxpair


def rev_compl(seq):
    """
    Does what you think it does.
    Handles the markup for junction type gracefully.
    """
    reversed_seq = seq[::-1]
    reversed_seq = reversed_seq.upper()
    rev_compl = ""
    for char in reversed_seq:
        if char == "A": rev_compl = rev_compl + "T"
        if char == "T": rev_compl = rev_compl + "A"
        if char == "G": rev_compl = rev_compl + "C"
        if char == "C": rev_compl = rev_compl + "G"
        if char == ":": rev_compl = rev_compl + ":"
        if char == ".": rev_compl = rev_compl + "."
        if char == "-": rev_compl = rev_compl + "-"
    return rev_compl


def consistent_maxpairs(list_of_maxpairs, bkpt, minlen=2):
    """
    takes a list of maxpairs, each in the form (index1, index2, wordlen, orientation)
    and returns a list containing the subset of maxpairs in the input list that are
    SD-MMEJ consistent and have repeated motifs of at least minlen bp.
    """
    consistent_maxpairs = []

    for maxpair in list_of_maxpairs:
        leftindex, rightindex, wordlen, ori = maxpair[0], maxpair[1], maxpair[2], maxpair[3]

        #forward repeats sitting directly on top of themselves
        #cannot be SD-MMEJ consistent
        
        if leftindex == rightindex and ori == "fwd":
            continue

        #rc repeats sitting on top of themselves can be SD-MMEJ consistent if
        #the breakpoint does NOT bisect the repeat (they're zero-loop hairpins)

        if leftindex == rightindex and ori == "rc":
            if bkpt == leftindex + (wordlen/2): continue
        
        #a repeat is SD-MMEJ consistent iff at least one
        #repeated motif includes the bkpt as an internal index
            
        if bkpt > leftindex and bkpt < (leftindex + wordlen):
            consistent_maxpairs.append(maxpair)
        elif bkpt > rightindex and bkpt < (rightindex + wordlen):
            consistent_maxpairs.append(maxpair)
            
    return consistent_maxpairs


def is_microhomology_jnxn(seq, deln, leftnick, rightnick):
    """
    seq is a string.
    deln is a tuple (L,R) of slice indices indicating the
    deletion boundary.
    is a microhomology junction if
    EITHER:
    *the rightmost base in the left-hand remaining chunk
    is the same as the rightmost base in the deleted chunk
    *the leftmost base in the right-hand remaining chunk is
    the same as the leftmost base in the deleted chunk
    AND
    * the ambiguous span is NOT completely to the left or
    completely to the right of both nicks.

    example - this is NOT a microhomology junction
    if the nicks are (0,4):
    TTATCCCTA
    TTAT--CTA #deln_indices = (4,6)

    this IS a microhomology junction:
    seq = "TTACCCCTA"
    leftnick, rightnick = 0,4
    deln_indices = (4,5)
    corresponds to
    TTACCCCTA
    TTAC-CCTA

    seq = "TTATCCCTA"
    deln = (4,6) #corresponds to TTAT--CTA
    leftnick, rightnick = 0, 4
    print is_microhomology_jnxn(seq, deln, leftnick, rightnick)
    #False

    seq = "TTACCCCTA"
    deln = (4,6) #corresponds to TTAC--CTA
    leftnick, rightnick = 0, 4
    print is_microhomology_jnxn(seq, deln, leftnick, rightnick)
    #True

    seq = "TTTTCCCTA"
    deln = (2,4) #corresponds to TT--CCCTA
    #this is a case of microhomology within overhangs
    leftnick, rightnick = 0, 4
    print is_microhomology_jnxn(seq, deln, leftnick, rightnick)
    #True

    seq = "TTTTCCCCAAAA"
    deln = (2,4) #corresponds to TT--CCCCAAAA
    leftnick, rightnick = 4, 8
    print is_microhomology_jnxn(seq, deln, leftnick, rightnick)
    #False
    
    """
    Ldel, Rdel = deln
    leftdna, deleted, rightdna = seq[:Ldel], seq[Ldel:Rdel], seq[Rdel:]
    mh_indices = id_jnxnal_mh(seq, deln, leftnick, rightnick)
    if mh_indices != None:
        ((LL, LR), (RL, RR)) = mh_indices
        #?rightmost in leftdna same as rightmost in deleted, and
        #?ambiguous span is not completely to the left of the left nick?
        if leftdna[len(leftdna)-1:] == deleted[len(deleted)-1:] and RR > leftnick:
            return True
        #?leftmost in rightdna same as leftmost in deleted, and
        #ambiguous span is not completely to the right of the right nick?
        if rightdna[0] == deleted[0] and LL < rightnick:
            return True
    else:
        return False


def get_jnxn_type(seq, deln, leftnick, rightnick):
    """
    helper function for EJR_junction constructor
    """
    if is_microhomology_jnxn(seq, deln, leftnick, rightnick) == False: return "ABJ"
    else: return "MHJ"


def push_all_mh_left(seq, deln_indices):
    """
    auxfn for id_jnxnal_mh
    WARNING: all this does is shift ambiguous bases flanking
    a deletion.  it does NOT check to make sure that the repeated
    motif spans a breakpoint - that is done within id_jnxnal_mh.
    seq = "CAATTTAAAAG"
    deln_indices = (1,7)
    corresponds to this alignment:
    CAATTTAAAAG
    C------AAAG (1,7)
    which should be
    CAATTTAAAAG
    CAA------AG (3,9)
    print push_all_mh_left(seq, deln_indices)
    (3, 9)
    """
    L,R = deln_indices
    while seq[L] == seq[R]:
        (L,R) = (L+1, R+1)
    return (L,R)


def id_jnxnal_mh(seq, deln_indices, leftnick, rightnick):
    """
    takes a sequence and a tuple of slice indices representing
    deletion boundaries and the left/right nick in the
    original sequence.  returns a tuple of the indices of the repeated   
    motifs in the original if there is junctional microhomology,
    None otherwise.
    """
    #convert deln indices to assign the ambiguous bases to the
    #left of the deletion:
    deln_indices = push_all_mh_left(seq, deln_indices)
    (L,R)=deln_indices
    left_mh, right_mh = [L,L], [R,R] 
    left_test, right_test = L-1, R-1
    leftdna, deleted, rightdna = seq[:L], seq[L:R], seq[R:]
    #check the base to the left of the left mh boundaries and move the
    #left mh boundaries left until the bases immediately to the left of
    #the left mh boundaries do not match
    while seq[left_test] == seq[right_test]:
        left_mh[0] -=1
        right_mh[0] -=1
        left_test -=1
        right_test -=1
    if left_mh[1] - left_mh[0] == 0:
        return None  #no junctional mh
    else:
        mh_indices = ((left_mh[0], left_mh[1]), (right_mh[0], right_mh[1]))
        ((LL, LR), (RL, RR)) = mh_indices
        if LL >= rightnick: return None
        elif RR <= leftnick: return None
        else: return mh_indices



def mh_len(mh_indices):
    """
    helper function
    mh_indices is the output of id_jnxnal_mh(seq, deln_indices)
    which is a tuple of tuples of slice indices of the microhomology
    repeats in the original sequence, ie
    mh_indices = ((1, 3), (7, 9))
    print mh_len(mh_indices)
    2
    """
    if mh_indices == None: return 0
    x,y = mh_indices[0][0], mh_indices[0][1]
    assert y>= x
    return y-x



def all_consistent_repeats(seq, deln, leftnick, rightnick):
    """
    returns a list of all SD-MMEJ consistent repeats for a given deletion
    """
    L, R = deln
    if not is_microhomology_jnxn(seq, deln, leftnick, rightnick):
        bkpt1 = L
        repaired_seq = seq[:L] + seq[R:]
        cm = consistent_maxpairs(unique_maxpairs(repaired_seq), bkpt1)
    if is_microhomology_jnxn(seq, deln, leftnick, rightnick):
        conv_deln_indices = push_all_mh_left(seq, deln)
        mh_indices = id_jnxnal_mh(seq, conv_deln_indices, leftnick, rightnick)
        bkpt1 = conv_deln_indices[0] - mh_len(mh_indices)
        bkpt2 = conv_deln_indices[0]
        repaired_seq = seq[:conv_deln_indices[0]] + seq[conv_deln_indices[1]:]
        cm = consistent_maxpairs(unique_maxpairs(repaired_seq), bkpt1)
        cm2 = consistent_maxpairs(unique_maxpairs(repaired_seq), bkpt2)
        for item in cm2:
            if item not in cm: cm.append(item) #don't double add reps that hit both ligation points
    return cm


def get_ligation_points(seq, jnxn_type, deln, leftnick, rightnick):
    """
    helper function for EJR_junction class:
    self.ligation_pointlist = get_ligation_points(self.original_seq, self.jnxn_type, raw_deln_indices)
    returns a list of 2 ligation points;
    for an ABJ, they will be the same - having the list be 2 items long
    for both junction types avoids having to special-case computations
    by junction type
    """
    ligation_points = []
    if jnxn_type == "ABJ":
        L, R = deln
        ligation_points.append(L)
        ligation_points.append(L)
        return ligation_points
    else:
        conv_deln_indices = push_all_mh_left(seq, deln)
        mh_indices = id_jnxnal_mh(seq, conv_deln_indices, leftnick, rightnick)
        bkpt1 = conv_deln_indices[0] - mh_len(mh_indices)
        bkpt2 = conv_deln_indices[0]
        ligation_points.append(bkpt1)
        ligation_points.append(bkpt2)
        return ligation_points


def get_motiflen(maxpair):
    """
    helper function for EJR_junction class
    returns the length of the repeated motif in a maxpair.
    reason this function is needed, instead of just using the
    wordlength, is because zero-loop hairpins are reported by the
    gstree that produces the maxpairs as one long sequence sitting on
    top of itself rather than as two adjacent motifs comprising an
    inverted repeat (zero-loop hairpins are their own reverse complements)
    """
    index1, index2, wordlen, ori = maxpair
    if index1 == index2: motiflen = wordlen/2
    else: motiflen = wordlen
    return motiflen


def get_deln_indices(ori_seq, deleted_seq):
    """
    takes an original sequence and a sequence with a single deleted
    interval indicated by dashes (-) and returns the slice indices of the
    deletion.  the original sequence must be given because this is
    meant as a helper function to extract deletion indices relative
    to the original sequence for the SD-MMEJ analysis jazz, and it needs
    to complain if the deletion is not indicated in the context of the
    original sequence, because then the slice indices will be wrong.
    the deletion can't run right up to either end of the deleted_seq;
    the first and last chars of deleted_seq must be a dna base.
    does NOT push all the microhomologies left.

    ori_seq =     "ATTACCCTGTTATCCCTAGCGGCCGCATAGGCC"
    deleted_seq = "ATTACCCTGTTATC----GCGGCCGCATAGGCC"
    print get_deln_indices(ori_seq, deleted_seq)
    #(14, 18)

    ori_seq =     "ATTACCCTGTTATCCCTAGCGGCCGCATAGGCC"
    deleted_seq = "ATTACCCTGTTATC-----GCGGCCGCATAGGCC"
    print get_deln_indices(ori_seq, deleted_seq)
    #original seq and deletion not properly aligned(A)

    ori_seq =     "ATTACCGTGTTATCCCTAGCGGCCGCATAGGCC"
    deleted_seq = "ATTACCCTGTTATC----GCGGCCGCATAGGCC"
    print get_deln_indices(ori_seq, deleted_seq)
    #original seq and deletion not properly aligned(B)

    ori_seq =     "ATTACCCTGTTATCCCTAGCGGCCGCTTAGGCC"
    deleted_seq = "ATTACCCTGTTATC----GCGGCCGCATAGGCC"
    print get_deln_indices(ori_seq, deleted_seq)
    #original seq and deletion not properly aligned(C)
    """
    if len(ori_seq) != len(deleted_seq):
        print "original seq and deletion not properly aligned(A)"
        return None
    ori_seq = ori_seq.upper()
    deleted_seq = deleted_seq.upper()
    j = 0
    while ori_seq[j] == deleted_seq[j]:
        j += 1
    if deleted_seq[j] == "-":
        del_index1 = j
    else:
        print "original seq and deletion not properly aligned(B)"
        return None
    while deleted_seq[j] == "-":
        j += 1
    if ori_seq[j] == deleted_seq[j]:
        del_index2 = j
    while ori_seq[j] == deleted_seq[j] and j < len(ori_seq)-1:
        j += 1
    if j < len(ori_seq)-1:
        print "original seq and deletion not properly aligned(C)"
        return None
    return (del_index1, del_index2)

    
        
        




