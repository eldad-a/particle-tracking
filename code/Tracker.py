#! python

################################################################################
#  filename: Tracker.py
#  
#  Copyright (c) 2014, Eldad Afik
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#  
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#  
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  
#  * Neither the name of this software nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
################################################################################

import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

# TODO:
#   * Need to verify that ghosts are handled correctly, namely - as future
#   prediction does not account for time difference, this is expected to result
#   in over-estimation of dX!

class Tracker():
    
    def __init__(self, positions_df, 
                 max_displacement=10, method=1, memory=0, 
                 condition_on_future=False, 
                 tracks=[], active_tracks=[], 
                 time_col='timestamp', position_cols=['x','y','z'], 
                 #start_frame=0, 
                 cleanup = True,
                 bounds=[], 
                 verbose=False):
        '''
        '''
        assert positions_df.index.name == time_col,\
         'expected positions_df to a have a single timestamps index' 
        self.positions = positions_df
        self.max_disp = max_displacement
        self.method = method # TODO: change the name (n_kinematic?)
        self.memory = memory
        self.condition_on_future = condition_on_future
        self.tracks = tracks # for each track holds [(timestamp, indx),pos]
        self.active_tracks = active_tracks
        assert time_col in positions_df.index.names,\
                'Could not find %s in the positions DataFrame index' % time_col
        self.time_col = time_col # must be a single int 
        self.pos_cols = [col for col in position_cols 
                             if col in positions_df.columns]
        #self.Ndim = self.positions[:, self.pos_cols].shape[1]
        self.Ndim = len(self.pos_cols) 
        assert self.Ndim>1, 'Found only %d position columns' % Ndim 
        self.frame_indx = 0 #start_frame   # not verified yet...!
        # if cleanup=True: upon analysing the last frame perform cleanup
        # assuming no more frames are to be analysed (such as removing ghosts
        # from the tracks)
        self.cleanup = cleanup 
        # TODO: verify that the treatment of the bounds is done correctly...
        self.bounds = bounds
        self.verbose = verbose
        if verbose:
            self.x = 0
            self.v = 0
            self.a = 0 
            print 'Time col is the index named %s' % self.time_col
            print '%d position columns were found, named: \n %s'\
                    % (self.Ndim, self.pos_cols)

    def user_update(self, matched, links):
        print '\nProcessed frame %d / %d ' % (self.frame_indx+1, self.timestamps.size)
        print '    Number of particles found:  %7d' % len(matched)
        print '    Number of active tracks:  %7d' % len(self.active_tracks)
        print '    Number of tracks with ghosts:  %7d' % len(self.ghost_tracks.keys())
        print '    Number of new tracks:  %7d' % np.sum(matched==0)
        print '    Number of terminated tracks:  %7d' % np.sum(links==-1)
        print '    Total number of tracks:  %7d' % len(self.tracks)
    
    def update_tracks(self, links, frame):
        """
        """
        t = self.timestamps[self.frame_indx]
        # attach the matched particles to their tracks
        matched = np.zeros(frame.shape[0], np.int64)
        for active_track_indx, track_indx in enumerate(self.active_tracks): 
            if links[active_track_indx]!=-1: # in case a link was found
                # if the ghost counter is non-zero, delete the ghosts from the
                # track
                try:
                    Nghosts = self.ghost_tracks[track_indx]
                    del self.tracks[track_indx][-Nghosts:]
                    del self.ghost_tracks[track_indx]
                except KeyError:
                    pass # no ghosts on this track      
                indx = links[active_track_indx]
                pos = frame[indx]
                self.tracks[track_indx].append(((t,indx),pos))
                matched[links[active_track_indx]]=1
            else:
                # If the ghost counter has not exceeded the memory parameter
                # append the esitmated position to the track, update the ghosts
                # counter and switch links to a non -1 such that the track
                # remains active (this is valid for non-singletons only,
                # there's no point in appending ghosts to singletons...)
                # if the ghost counter exceeds the memory value, remove the
                # ghosts from the track

                ## TODO: once estimation will be based on timestamp as well,
                ## will need to convert it to a pd.Series with a timestamp
                ghost = (None,self.estimate[active_track_indx])
                try:
                    if self.ghost_tracks[track_indx] < self.memory:
                        self.tracks[track_indx].append(ghost)
                        self.ghost_tracks[track_indx] += 1
                        links[active_track_indx] = -2
                    else:
                        Nghosts = self.ghost_tracks[track_indx]
                        del self.tracks[track_indx][-Nghosts:]
                        del self.ghost_tracks[track_indx]
                except KeyError:
                    if self.memory and len(self.tracks[track_indx])>1:
                        self.tracks[track_indx].append(ghost)
                        self.ghost_tracks[track_indx] = 1
                        links[active_track_indx] = -2
                        
        # update active_tracks (terminate those which did not link)
        self.active_tracks = [
                        track_indx for active_track_indx, track_indx in 
                             enumerate(self.active_tracks) 
                                    if links[active_track_indx]!=-1
                             ]
        # start new tracks for the unmatched particles
        first = len(self.tracks)
        self.active_tracks += range(first, first+sum(matched==0))
        unmatched = np.nonzero(matched==0)[0]
        self.tracks += [    
                        [((t,unmatched[indx]),pos)] 
                             for indx,pos in enumerate(frame[unmatched])
                       ]
        if self.verbose: self.user_update(matched, links)
        
    def estimate_future_position(self, Xt):
        """
        Estimate future position based on (possibly) three known past positions
        by parabolic extrapolation (equivalent to Taylor expansion neglecting
        order dt^3).
        If there are not enough past points attempt linear approximation (2
        points).
        For singletons use current position (similar to diffusion).
        """
        pos = 1
        # TODO: account for the inter-frame time difference when available
        try:
            ## The following can be found by either Taylor expansion up to dt^3
            ## or equivalently by parabolic extrapolation
            estimate = 3*Xt[-1][pos] - 3*Xt[-2][pos] + Xt[-3][pos]
            # Seems like N.Ouellette is using a weighted parabolic
            # extrapolation:
            #estimate = 2.5*Xt[-1][pos_cols] \
            #         - 2*Xt[-2][pos_cols] \
            #         + .5*Xt[-3][pos_cols]
            if self.verbose: self.a += 1
            return estimate
        except IndexError: # trajectory contains less than 3 points:
            pass # fallback to linear extrapolation
        try:
            ## Use linear extrapolation:
            estimate = 2*Xt[-1][pos] - Xt[-2][pos]
            if self.verbose: self.v += 1
        except IndexError: # trajectory contains less than 2 points:
            estimate = Xt[-1][pos] # use the current position 
            if self.verbose: self.x += 1
        return estimate

    def get_cost_mtx(self, next_pos):
        """
        """
        # TODO: account for the inter-frame time difference when available
        self.estimate = np.empty((len(self.active_tracks), self.Ndim))
        for m,track_indx in enumerate(self.active_tracks): 
            # TODO: can this loop be turned into an array op?
            track = self.tracks[track_indx]
            self.estimate[m] = \
                    self.estimate_future_position(track[-1-self.method:])
        #try:
        #    return cdist(self.estimate, next_pos, 'sqeuclidean')
        #except ValueError as err:
        #    # assuming err.message == 'XB must be a 2-dimensional array.'
        #    next_pos = next_pos[np.newaxis] ## convert to a 2d array
        #    return cdist(self.estimate, next_pos, 'sqeuclidean')
        return cdist(self.estimate, next_pos, 'sqeuclidean')
    
    def one_frame_linking(self):
        """
        Simple one frame linking based on best match between estimated position
        of tracer to the measured positions of tracers in the next frame.
        Best match is defined by smallest (Euclidean) distance.
        Future position estimate is based on extrapolation up to second order
        (three pre-linked positions) in inter-frame time difference.
        """
        # TODO: currently not checking/warning for empty frames
        t = self.timestamps[self.frame_indx]
        frame1 = self.positions.loc[t, self.pos_cols].values
        if frame1.ndim == 1:
            frame1 = frame1[np.newaxis]
        dist_mtx = self.get_cost_mtx(frame1)
        
        # initialise costs and links
        costs = np.zeros(len(self.active_tracks))
        links = -np.ones_like(costs, np.int64)
        
        # assuming the new frame contains particles
        for active_track_indx in xrange(len(self.active_tracks)): 
            # TODO: can this loop be turned into an array op?
            distances = dist_mtx[active_track_indx]
            costs[active_track_indx] = distances.min()
            ## end track if best estimate is out of the search radius:
            if costs[active_track_indx] > self.max_disp2:
                continue
            ## are there two best estimates? in case of ambiguity end track:
            if np.sum(distances==costs[active_track_indx]) > 1: continue
            # check whether another track was already linked this particle
            bestmatch = distances.argmin()
            previous_match = links == bestmatch
            if previous_match.sum():
                if costs[previous_match] > costs[active_track_indx]: 
                    # if the link is "cheaper" undo previous
                    links[previous_match] = -1
                else:
                    continue # really?
            # found link!
            links[active_track_indx] = bestmatch
        self.update_tracks(links, frame1)

    def linker(self):
        '''
        expects positions as a pandas DataFrame, of position cols named as
        stored in pos_cols and a timestamp index with name as stored in
        time_col.
        if provided with tracks, continues appending to the already 
        exisiting trajectories, starting from the begining of positions.
        '''
        ###            Setting-up                     ###    
        self.max_disp2 = self.max_disp*self.max_disp
        self.timestamps = self.positions.index.\
                                            get_level_values(self.time_col).\
                                            unique()
        assert len(self.timestamps)>1, \
                'provided with one frame only, linking is trivial...'
        # In case memory is non zero, tracks which were not linked 
        # will be appended with estimated position ("ghost"), until
        # either they will link or the memory limit is reached.
        # This is done using a dictionary to hold the number of ghosts appended 
        # with the key holding the track_indx
        self.ghost_tracks = {} 
        
        if len(self.tracks): 
            assert False,\
                    'linking with history (pre-tracks) has not been tested yet'
        else:
            # in case no exisiting tracks were given need to initialise
            t = self.timestamps[0]
            frame0 = self.positions.loc[t,self.pos_cols].values
            if frame0.ndim == 1: 
                frame0 = frame0[np.newaxis]
            self.tracks = [
              #[((t,indx),pos)] for indx,(t,pos) in enumerate(frame0.iterrows())
              [((t,indx),pos)] for indx, pos in enumerate(frame0)
                          ]
            # activate all tracks (add all to active_tracks)
            self.active_tracks = range(len(self.tracks))
            matched = np.zeros(len(self.active_tracks))
            links = [0]
            if self.verbose: self.user_update(matched, links)
            self.frame_indx += 1
        
        ###            Linking main loop              ###
        for self.frame_indx in xrange(self.frame_indx, self.timestamps.size):
            if self.condition_on_future:
                self.two_frames_linking() # NO & HX 4BE
            else:
                self.one_frame_linking()
        
        # delete remaining ghosts from tracks:
        if self.cleanup:
            for track_indx, Nghosts in self.ghost_tracks.iteritems():
                del self.tracks[track_indx][-Nghosts:]
                #del self.ghost_tracks[track_indx]             
            
        # TODO prune short tracks?
        # TODO calculate tracks statistics (such as lengths)
        # TODO visualisation

    def update_traj_ids(self, min_length=2, traj_id=0, unlinked_id=0):
        '''
        Tracers in tracks shorter than min_length will be discarded (default 2).
        Tracer which do not belong to any track are assigned with unlinked_id
        (default 0). 
        The first traj_id will be the input traj_id + 1 (default traj_id=0)
        '''
        ## post-process tracks - get a DataFrame for each and a list of indices
        ## of the tracers index in the frame (in the input positions_df)
        summary = []
        self.last_traj_id = traj_id
        ## get grouped indices frames:
        ts_indices = self.positions.groupby(level=self.time_col).indices
        ## for better performance, in updating avoid the DataFrame indexing.
        # (do this buy getting a view of the `traj_id` column and update by
        # numpy indexing)
        try:
            traj_ids = self.positions['traj_id'].values
        except KeyError:
            #traj_ids = np.empty_like(self.positions[self.positions.columns[0]])
            #traj_ids = unlinked_id
            self.positions['traj_id'] = unlinked_id
            traj_ids = self.positions['traj_id'].values
            if self.verbose:
                print 'Initialised as "traj_id" column to %s' % unlinked_id
        ## update table:
        for track in self.tracks:
            if len(track)<min_length: continue # skip short tracks 
            traj_id+=1
            track_indices, track = zip(*track)
            to_update = [ts_indices[t][indx] for (t,indx) in track_indices]
            #self.positions.ix[to_update,'traj_id'] = traj_id ## slow!
            #self.positions['traj_id'].iloc[to_update] = traj_id ## faster
            traj_ids[to_update] = traj_id # update directly on the numpy view
            ti,tf = track_indices[0][0], track_indices[-1][0]
            Nframes = len(track_indices)
            summary.append((traj_id, ti, tf, Nframes))
            self.last_traj_id = traj_id
        summary = pd.DataFrame.from_records(summary,
                                  columns=['traj_id','Ti','Tf','Nframes'], 
                                  index='traj_id'
                                            )
        self.summary = summary
        if self.verbose:
            print "Found %d tracks longer than %d." %\
                                                (len(summary),min_length)


################# experimental
'''
DISCLAIMER:
This part should be regarded as incomplete and highly experimental.
It is provided ONLY in order  to help others who may like to implement the algorithms described in

N. T. Ouellette, H. Xu, and E. Bodenschatz, 
A quantitative study of three-dimensional Lagrangian particle tracking algorithms,
Exp. Fluids 40, 301-313 (2006).

My personal experience shows that the output based on this approach is not as good as the ones by the main one.
Use with care!
'''

    def two_frames_linking(self):
        """
        """
        err =  'conditioning on future linking is not implemented yet '
        err += 'for pandas exiting...'

        assert False, err
        # TODO: currently not checking/warning for empty frames
        
        ## Load the next and next-next future frames
        # TODO: This way all (but boundary) frames are to be loaded twice.
        # need to verify that this is not wasting time (as this could be done differently)
        frame1 = self.positions.xs(self.timestamps[self.frame_indx],
                                   #level=self.time_col, 
                                   copy=False,
                                   )
        try:
            frame2 = self.positions.xs(self.timestamps[self.frame_indx+1],
                                   #level=self.time_col, 
                                   copy=False,
                                   )
        except IndexError:
            # in case this the next frame is the last one revert to single frame linking:
            self.one_frame_linking()
            return
        
        # initialise costs and links
        costs = np.zeros(len(self.active_tracks))
        links = -np.ones_like(costs)
        
        # get all potential links within search radius
        # for all active_tracks:
        #    estimate future position based on potential links
        #    choose the potential link which results in smallest distance from next future positions
        dist_mtx = self.get_cost_mtx(frame1)
        #potential_links = np.transpose(np.nonzero(dist_mtx <= self.max_disp2))
        # TODO: account for the inter-frame time difference when available
        for active_track_indx, distances in enumerate(dist_mtx):
            past_positions = self.tracks[self.active_tracks[active_track_indx]][-self.method:]
            potential_links = distances <= self.max_disp2
            if not potential_links.sum():
                # in case there are no potential links for this track, move on to the next track
                costs[active_track_indx] = distances.min() # put the negative of the minimal distance for debug
                continue
            estimate = np.empty((potential_links.sum(), self.Ndim))
            for m,potential_position in enumerate(frame1[potential_links]):
                # use "+" operator for list concatenation (TODO: be cautious!)
                estimate[m] = self.estimate_future_position(past_positions+[potential_position])
            potential_costs = cdist(estimate, frame2[:, self.pos_cols], 'sqeuclidean')
            costs[active_track_indx] = potential_costs.min()
            ## end track if best estimate is out of the search radius:
            if costs[active_track_indx] > self.max_disp2: continue
            potential_links_indices = np.nonzero(potential_links)[0]
            minimal_cost_potential_links = np.nonzero(potential_costs==costs[active_track_indx])[0]
            bestmatch_indx = potential_links_indices[minimal_cost_potential_links]
            # TODO: consider avoiding the following step
            if len(bestmatch_indx) > 1: # are there two best estimates? 
                print '\n resolving ambiguity!'
                ## in case of ambiguity consider also the earlier costs and choose the minimal among them:
                bestmatch_indx = bestmatch_indx[distances[bestmatch_indx].argmin()]
                # are there two best estimates? in case of ambiguity consider terminating track
                # (currently choose the first one) 
            ## resolve multiple assignments to the same particle:
            previous_match = links == bestmatch_indx
            if previous_match.sum():
                if costs[previous_match] > costs[active_track_indx]: # indeed bestmatch
                    # TODO: attempt simple bestmatch
                    links[previous_match] = -1
                else:
                    # TODO: attempt simple bestmatch
                    continue # really?
            links[active_track_indx] = bestmatch_indx
            #print 'linked %s -> %s' %\
            #       (self.tracks[self.active_tracks[active_track_indx]][-1], frame1[links[active_track_indx]])  
        self.update_tracks(links, frame1)

    def __one_frame_linking_two_frames_style__(self):
        """
        """
        err =  'conditioning on future linking is not implemented yet '
        err += 'for pandas exiting...'
        assert False, err
        # TODO: currently not checking/warning for empty frames
        
        ## Load the next and next-next future frames
        # TODO: This way all (but boundary) frames are to be loaded twice.
        # need to verify that this is not wasting time (as this could be done differently)
        frame1 = self.positions[self.frames_indices[self.frame_indx]:\
                                self.frames_indices[self.frame_indx+1]]
        # initialise costs and links
        costs = np.zeros(len(self.active_tracks))
        links = -np.ones_like(costs)
        
        # get all potential links within search radius
        # for all active_tracks:
        #    estimate future position based on potential links
        #    choose the potential link which results in smallest distance from next future positions
        dist_mtx = self.get_cost_mtx(frame1)
        #potential_links = np.transpose(np.nonzero(dist_mtx <= self.max_disp2))
        # TODO: account for the inter-frame time difference when available
        
        for active_track_indx, distances in enumerate(dist_mtx):
            past_positions = self.tracks[self.active_tracks[active_track_indx]][-1-self.method:]
            potential_links = distances <= self.max_disp2
            if not potential_links.sum():
                # in case there are no potential links for this track, move on to the next track
                costs[active_track_indx] = distances.min() # put the negative of the minimal distance for debug
                continue
            estimate = np.empty((1, self.Ndim))#np.empty((potential_links.sum(), self.Ndim))
            if True:#for m,potential_position in enumerate(frame1[potential_links]):
                m=0
                # use "+" operator for list concatenation
                potential_position = past_positions[-1]
                tst = [potential_position] if len(past_positions)==1 else past_positions[:-1]+[potential_position]
                estimate[m] = self.estimate_future_position(tst)
            potential_costs = cdist(estimate, frame1[:, self.pos_cols], 'sqeuclidean')
            costs[active_track_indx] = potential_costs.min()
            ## end track if best estimate is out of the search radius:
            if costs[active_track_indx] > self.max_disp2: continue
            potential_links_indices = np.nonzero(potential_links)[0]
            minimal_cost_potential_links = np.nonzero(potential_costs==costs[active_track_indx])[1]
            bestmatch_indx = minimal_cost_potential_links#potential_links_indices[minimal_cost_potential_links]
            # TODO: consider avoiding the following step
            if len(bestmatch_indx) > 1: # are there two best estimates? 
                print '\n resolving ambiguity!'
                ## in case of ambiguity consider also the earlier costs and choose the minimal among them:
                bestmatch_indx = bestmatch_indx[distances[bestmatch_indx].argmin()]
                # are there two best estimates? in case of ambiguity consider terminating track
                # (currently choose the first one) 
            ## resolve multiple assignments to the same particle:
            previous_match = links == bestmatch_indx
            if previous_match.sum():
                if costs[previous_match] > costs[active_track_indx]: # indeed bestmatch
                    # TODO: attempt simple bestmatch
                    links[previous_match] = -1
                else:
                    # TODO: attempt simple bestmatch
                    continue # really?
            links[active_track_indx] = bestmatch_indx
            #print 'linked %s -> %s' %\
            #       (self.tracks[self.active_tracks[active_track_indx]][-1], frame1[links[active_track_indx]])
        self.update_tracks(links, frame1)
