#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 20:39:12 2020

@author: transcend
"""


import re
import numpy as np
import pandas as pd
import random


class robanalysis(object):
    def __init__(self):
        
        """Object attributes
        self.Reg_ss: Vector of the single site regulatory sequences (15bp)
        
        self.Reg_os: Vector of the overlapping site regulatory sequences (15bp)
        
        self.PWM_Kdref: PWM in the format of Kd/Kdref, so that to obtain the Kd 
        of a given sequence, you can multiply all the Kd/Kdref factors based on 
        positions in the PWM, and then multiply this by the Kd of the consensus
        
        self.Kd_consensus: Kd of the consensus sequence
        
        self.ss_affinity_df: dataframe for single site, with row headers corresponding
        to the number of mutations, and columns corresponding to the trial# (starting
        at zero... so if 50 trials, then labelled from 0 to 49)
        
        self.os_affinity_df: overlapping site equivalent
        """
        
        self.Nr_trials = None
        self.Reg_ss = None
        self.Reg_os = None 
        self.PWM_Kdref = None
        self.Kd_consensus = None
        self.ss_affinity_df = None
        self.os_affinity_df = None
        
    def launch_simulation(self, Kd_cons, Nr_mutations):
        """
        This is the main function that iterates through:
            >PWM scanning the regulatory sequences by calling __pwm_scan, 
            >multiplying the affinity vectors output by the Kd of the consensus
            >storing the affinity vector values as columns of an affinity dataframe
            >calling a function to randomly mutate a single base of each regulatory sequence 
        
        Args:
            Kd_consensus: Kd of the consensus sequence.
            Nr_mutations: corresponds to the number of iterations to run through

        Returns:
            fills out the self.ss_affinity_df and self.os_affinity_df:
        """
        
        self.Kd_consensus = Kd_cons
             
        self.ss_affinity_df = pd.DataFrame(columns=range(0,self.Nr_trials)) #create dataframe to store the results of simulation (rows = #mutations, columns=sequence#)
        
        self.os_affinity_df = pd.DataFrame(columns=range(0,self.Nr_trials))
        
        self.ss_affinity_df = self.ss_affinity_df.append(pd.Series(Kd_cons*self.__pwm_scan(self.Reg_ss)),ignore_index=True) #Store the intial affinity values (should all be the same), before any mutation
                                                                                                       #note that these columns will have row-index 0 (signifying 0 mutations)
        self.os_affinity_df = self.os_affinity_df.append(pd.Series(Kd_cons*self.__pwm_scan(self.Reg_os)),ignore_index=True)
        
        for i in range(1,Nr_mutations+1):
            
            #randomly substitution-mutate each sequence in the sequence vectors
            self.Reg_ss = self.substitution_mutate(self.Reg_ss) #seems instinctive to copy the vector, mutate the copy, and then update the original using the mutated copy...
            
            self.Reg_os = self.substitution_mutate(self.Reg_os)
            
            #check the affinities of the new sequence vectors, and record them in the dataframe
            self.ss_affinity_df = self.ss_affinity_df.append(pd.Series(Kd_cons*self.__pwm_scan(self.Reg_ss)),ignore_index=True)
            
            self.os_affinity_df = self.os_affinity_df.append(pd.Series(Kd_cons*self.__pwm_scan(self.Reg_os)),ignore_index=True)
            
            
        
    def load_kdref_pwm(self, filename, n_mer):
        """
        Args: 
            filename: corresponds to a file that contains a Kd/Kdref PWM (ex: from Swank 2019
            num_bases: number of bases in the motif
        """
        PWM_Kdref = np.ones((4, n_mer), np.float) #Start with 1s, because Kdref/Kdref = 1
        
        with open(filename, 'r') as fh:
            PWM_Kdref = [[float(B) for B in line.split()] for line in fh] #process the text file into a 2D list
        
        PWM_Kdref = np.array(PWM_Kdref, dtype=float) #convert the 2D array into a 2D numpy array
        
        PWM_Kdref[PWM_Kdref < 1] = 1 #minimum value is 1
                
        self.PWM_Kdref = PWM_Kdref
        
        return
        
    def load_sequence(self, seq_ss, seq_os, nr_trials):
        """
        Args:
            seq_ss: the single-site DNA sequence, containing only (A, C, G, T)
                    no space or other characters allowed
            seq_os: the overlapping-site DNA sequence
            nr_trials: number of trials of regulatory sequences to analyze
        """
        self.Nr_trials = nr_trials
        
        numeric_ss = self.__str_to_np_seq(seq_ss) # Convert the single site sequence into numeric encoding (a np.array)
        numeric_os = self.__str_to_np_seq(seq_os) # Convert the overlapping site sequence into numeric encoding (a np.array)
        
        self.Reg_ss = np.tile(numeric_ss, (nr_trials, 1)) # Create vector of the single site sequence
        self.Reg_os = np.tile(numeric_os, (nr_trials, 1)) # Create vector of the overlapping site sequence
        
        return

    def __str_to_np_seq(self, str_seq):
        """
        A custom DNA base coding system with numbers.
        (A, C, G, T, N) = (0, 1, 2, 3, 0)

        A DNA string is converted to a numpy integer array (np.unit8) with the same length.
        """
        np_seq = np.zeros((len(str_seq), ), np.uint8)

        ref = {'A':0, 'a':0,
               'C':1, 'c':1,
               'G':2, 'g':2,
               'T':3, 't':3,}

        for i, base in enumerate(str_seq): #i is the index, base is the letter  in str_seqs
            np_seq[i] = ref.get(base) #ref is a dictionary, using get allows to use base as a key to get associated numeric value

        return np_seq

    def __pwm_scan(self, seq_vec):
        """
        The core function that performs PWM scan
        through the sequence and returns a vector of the values of the highest affinity sites.
        """

        n_mer = self.PWM_Kdref.shape[1]    # length (num of cols) of the weight matrix
        
        cols = np.arange(n_mer) # column indices for PWM from 0 to (n_mer-1)
        
        score = 50000
        
        PWM_Kdref_rc = self.PWM_Kdref[::-1, ::-1] # Reverse complementary PWM
        
        vector_affinities = np.zeros((seq_vec.shape[0]),float)
        
        
        # The main loop that scans through the regulatory sequences
        
        for j in range(seq_vec.shape[0]): #iterate over every member of Reg_ss and Reg_os
            
            for i in range(seq_vec.shape[1] - n_mer + 1): #iterate over every base in a given member of Reg_ss and Reg_os
    
                window = seq_vec[j,i:(i+n_mer)] #pull out the sequence we're comparing the PWM against.
    
                # --- The most important line of code ---
                #     Use integer coding to index the correct score from column 0 to (n_mer-1)
                new_score = np.prod( self.PWM_Kdref[window, cols] ) #indexing with two np arrays creates a new array where members are 
                                                    #specified as ([rows: (first member's row index), (second member's row index), etc.],
                                                    #[columns: (first member's col index), (second member's col index), etc.])
                                                    #note that the score here is the Kd/Kdref, we will preserve the highest score
                                                    #and then before registering it in the dataframe, multiply it by the Kd
                if new_score < score:   #keep only the highest score 
                    score = new_score
                
                new_score = np.prod( PWM_Kdref_rc[window, cols] ) #technically it should not matter if the highest affinity site
                                                                       #is on the other strand --> only one site will dominate in this model
                                                                       #note: since seeding pre-sites on forward strand, unlikely other one will
                                                                        #win, but this is more general.
                if new_score < score:
                    score = new_score
                    
            vector_affinities[j] = score #record the lowest Kd/Kdref for the jth regulatory sequence 
            
            score = 50000 #reset the score to something unattainable
            
        return vector_affinities
    
    
    def substitution_mutate(self, seq_vec):
        
        copy_seq_vec = np.copy(seq_vec)
        
        for j in range(copy_seq_vec.shape[0]):
            
            rand_ind = random.randint(0,copy_seq_vec.shape[1]-1)
            
            if copy_seq_vec[j,rand_ind] == 0:
                copy_seq_vec[j,rand_ind] = random.choice([1,2,3])
                
            elif copy_seq_vec[j,rand_ind] == 1:
                copy_seq_vec[j,rand_ind] = random.choice([0,2,3])
                
            elif copy_seq_vec[j,rand_ind] == 2:
                copy_seq_vec[j,rand_ind] = random.choice([0,1,3])
                
            elif copy_seq_vec[j,rand_ind] == 3:
                copy_seq_vec[j,rand_ind] = random.choice([0,1,2])
                
        return copy_seq_vec
        
        
        
        
        