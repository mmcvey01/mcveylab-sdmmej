from EJR_junction_auxfns import *

class EJR_junction(object):
    """
    A container for generating and gluing together
    and a sieve for filtering information
    about an end-joining repair junction.
    """

#plan is to run all_valid_delns first and then use its output to create
#an EJR_junction object for each deletion.

    def __init__(self, original_seq, raw_deln_indices, leftnick, rightnick):
        """
        """
        self.raw_deln_indices = raw_deln_indices
        self.original_seq = original_seq
        self.repaired_seq = original_seq[:raw_deln_indices[0]] + original_seq[raw_deln_indices[1]:]
        self.jnxn_type = get_jnxn_type(original_seq, raw_deln_indices, leftnick, rightnick)
        self.converted_deln_indices = push_all_mh_left(original_seq, raw_deln_indices)
        #NOTE: self.converted_deln_indices depends on other code to make sure that this only
        #ever gets applied to microhomology junctions
        self.deln_size = raw_deln_indices[1] - raw_deln_indices[0]
        self.mh_rep_in_original_seq = id_jnxnal_mh(original_seq, raw_deln_indices, leftnick, rightnick)
        self.mh_len = mh_len(id_jnxnal_mh(original_seq, raw_deln_indices, leftnick, rightnick))
        self.ligation_pointlist = get_ligation_points(self.original_seq, self.jnxn_type, raw_deln_indices, leftnick, rightnick)
        self.consistent_repeats = all_consistent_repeats(original_seq, raw_deln_indices, leftnick, rightnick)
        self.leftnick = leftnick
        self.rightnick = rightnick

    def __str__(self):
        """
        """
        return str(((self.raw_deln_indices), (self.original_seq)))



    def filtered_cr(self, minlen, search_radius):
        """
        filtered consistent repeats
        returns a list containing the subset of consistent repeats
        with wordlength >= minlen and with at least minlen bp of both
        halves of the repeat.
        uses wordlen/2 as minlen for zero bp hairpins
        """
        ###NEEDS TESTING###
        filtered_by_minlen = []
        for rep in self.consistent_repeats:
            index1, index2, wordlen, ori = rep
            if ori == "fwd" and wordlen >= minlen: filtered_by_minlen.append(rep)
            if ori == "rc":
                if index1 != index2 and wordlen >= minlen: filtered_by_minlen.append(rep)
                if index1 == index2:
                    if wordlen/2 >= minlen: filtered_by_minlen.append(rep)
        filtered_by_search_radius = []
        #need to test only for failure conditions
        for rep in filtered_by_minlen:
            #convert indices for zero nucleotide hairpins
            index1, index2, wordlen, ori = rep
            if index1 == index2:
                index2 = index1 + (wordlen/2)
                wordlen = wordlen/2
            #fail if result of subtracting the slice index of the left search
            #radius bound from the slice index that is the right bound of the
            #left-hand motif (ie index1 + wordlen) is less than minlen. 
            #index of left search radius bound is the slice index of the
            #left ligation point (self.ligation_pointlist[0]) - search_radius
            left_search_radius_index = self.ligation_pointlist[0] - search_radius
            if (index1 + wordlen) - left_search_radius_index < minlen:
                #print rep, "filtered: left motif outside  left search radius = ", left_search_radius_index
                continue
            # fail if result of subtracting the index of the left bound of the
            #right-hand repeated motif (ie, index2) from the index of the
            # right-hand search radius boundary is less than minlen
            right_search_radius_index = self.ligation_pointlist[1] + search_radius
            if right_search_radius_index - index2 < minlen:
                #print rep, "filtered: right motif outside right search radius = ", right_search_radius_index
                continue
            else:
                #print "C"
                filtered_by_search_radius.append(rep)
        return filtered_by_search_radius
 


    def align_deln_conv_indices(self):
        """
        returns the internally deleted sequence associated with the deletion object
        represented with dashes indicating deleted base pairs for easy alignment
        with the original sequence.  uses converted deln indices.
        """
        showdeln = ""
        for index, value in enumerate(self.original_seq):
            if index < self.converted_deln_indices[0]: showdeln = showdeln + value
            elif index >= self.converted_deln_indices[1]: showdeln = showdeln + value
            else: showdeln = showdeln + "-"
        return showdeln


    def align_deln_raw_indices(self):
        """
        returns the internally deleted sequence associated with the deletion object
        represented with dashes indicating deleted base pairs for easy alignment
        with the original sequence.  uses the raw deletion boundaries (ie ambiguous
        bases have not been shoved over to the left)
        """
        showdeln = ""
        for index, value in enumerate(self.original_seq):
            if index < self.raw_deln_indices[0]: showdeln = showdeln + value
            elif index >= self.raw_deln_indices[1]: showdeln = showdeln + value
            else: showdeln = showdeln + "-"
        return showdeln
    

    def print_jnxn(self):
        """
        prints the repair junction as a string with apparent points of
        ligation indicated by upper/lowercase transitions.
        """
        to_print = ""
        if self.jnxn_type == "ABJ":
            for index, value in enumerate(self.repaired_seq):
                if index < self.ligation_pointlist[0]: to_print = to_print + value
                else: to_print = to_print + value.lower()
        if self.jnxn_type == "MHJ":
             for index, value in enumerate(self.repaired_seq):
                 if index < self.ligation_pointlist[0]: to_print = to_print + value
                 elif index >= self.ligation_pointlist[1]: to_print = to_print + value
                 else: to_print = to_print + value.lower()
        return to_print



    def print_fcr_aligns(self, minlen, search_radius):
        """
        print filtered consistent repeat alignments
        """
        fcr = self.filtered_cr(minlen, search_radius)
        print "="*60
        print "raw deletion indices are", self.raw_deln_indices
        print "search radius length is", search_radius
        print "search interval indices are", self.ligation_pointlist[0] - search_radius, self.ligation_pointlist[1] + search_radius
        print "minimum consistent repeat length is", minlen
        print "complete consistent repeat list is", self.consistent_repeats
        print "filtered consistent repeat list is", fcr
        jnxn_seq = self.print_jnxn()
        if fcr == []:
            print
            print "no filtered consistent repeats for this deln"
            print self.original_seq
            print self.align_deln(), self.raw_deln_indices
            print jnxn_seq, self.converted_deln_indices
            print
        for maxpair in fcr:
            i1, i2, wordlen, ori = maxpair
            motif1 = ""
            for index, value in enumerate(jnxn_seq):
                if index in range(i1, i1+wordlen):
                    motif1 = motif1 + value
                else: motif1 = motif1 + "-"
            motif2=""
            for index, value in enumerate(jnxn_seq):
                if index in range(i2, i2+wordlen):
                    motif2 = motif2 + value
                else: motif2 = motif2 + "-"
            print
            print self.original_seq
            print self.align_deln(), self.raw_deln_indices
            print jnxn_seq, self.converted_deln_indices
            print motif1, maxpair, "motif length is", get_motiflen(maxpair)
            print motif2
            print
        
                
            
            
                
            
        
            
            
        
        
