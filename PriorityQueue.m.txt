classdef PriorityQueue
%------------------------------------------------------------------------------
% PriorityQueue
%
% Methods:
%   [pq]     = PriorityQueue( capacity )
%   [flag]   = isempty( pq )
%   [n]      = nactive( pq )
%   [pq]     = insert( pq, e )
%   [pq, e]  = pop( pq )
%   [pq]     = cancel( pq, id )
%   [pq_new] = purge( pq )
%   [pq_new] = expand( pq )
%              dump( pq, fid )
%              dumpactive( pq, fid )
%
% Notes:
% o This code is part of Homework 10, CE4121, Spring 2013.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Written by:
%   Dr. Randal J. Barnes
%   Department of Civil Engineering
%   University of Minnesota
%
% Version:
%   7 April 2013
%-------------------------------------------------------------------------------

    properties (SetAccess = private)
        top        = nan;           % index of the next (top) event
        count      = 0;             % number of events in queue (includes zombies)
        capacity   = 0;             % current storage capacity of queue.
        next       = nan(0,1);      % array of next events (linked list)
        eventArray = cell(0,1);     % cell array of events
    end
    
    methods
        %-----------------------------------------------------------------------
        % Constructor with an initial size of <capacity>.
        %-----------------------------------------------------------------------        
        function [pq] = PriorityQueue( capacity )
            pq.top        = nan;
            pq.count      = 0;
            pq.capacity   = capacity;
            pq.next       = nan(capacity,1);
            pq.eventArray = cell(capacity,1);
        end

        %-----------------------------------------------------------------------
        % Returns true if there are no active events in the priority queue.
        %-----------------------------------------------------------------------        
        function [flag] = isempty( pq )
            flag = isnan( pq.top );
        end

        %-----------------------------------------------------------------------
        % Returns the number of active events in the priority queue.
        %-----------------------------------------------------------------------        
        function [n] = nactive( pq )
            n = 0;
            j = pq.top;
                
            while isfinite(j)
                n = n + 1;
                j = pq.next(j);
            end
        end
        
        %-----------------------------------------------------------------------
        % Insert a new event in the priority queue.  The routine purges and
        % expands the allocated storage if needed.
        %-----------------------------------------------------------------------        
        function [pq] = insert( pq, e )
            % Make certain that the priority queue has space.
            if pq.count == pq.capacity
                pq = purge(pq);
                if pq.count == pq.capacity
                    pq = expand(pq);
                end
            end
                
            % Put the new event in the next open space.
            pq.count = pq.count + 1;
            pq.eventArray{pq.count} = e;
            
            % Determine the place in the priority queue.
            if isnan(pq.top)
                pq.top = pq.count;
            elseif isbefore( e, pq.eventArray{pq.top} )
                pq.next(pq.count) = pq.top;
                pq.top = pq.count;
            else
                i = pq.top;
                j = pq.next(i);
                
                while isfinite( j )
                    if isbefore( e, pq.eventArray{j} )
                        break;
                    else
                        i = j;
                        j = pq.next(j);
                    end
                end
                pq.next(i) = pq.count;                
                pq.next(pq.count) = j;
            end
        end

        %-----------------------------------------------------------------------
        % Pop the top event off the priority queue.
        %-----------------------------------------------------------------------        
        function [pq, e] = pop( pq )
            if isnan( pq.top )
                e = nan;
            else
                e = pq.eventArray{pq.top};
                pq.top = pq.next(pq.top);
            end
        end
 
        %-----------------------------------------------------------------------
        % Cancel all events involving <id>.  The cancelled events are not
        % erased from the priority queue, they are simply disconnected.
        %-----------------------------------------------------------------------
        function [pq] = cancel( pq, id )
            % Cancel the leading events that involve id.
            while isfinite( pq.top )
                if involves( pq.eventArray{pq.top}, id )
                    pq.top = pq.next( pq.top );
                else
                    break;
                end
            end

            % Cancel any events that involve bid.            
            if isfinite( pq.top )
                i = pq.top;
                while isfinite( pq.next(i) )
                    j = pq.next(i);
                    if involves( pq.eventArray{j}, id )
                        pq.next(i) = pq.next(j);
                    else
                        i = j;
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------------
        % Eliminate all zombie events and move all active events to the top
        % of the allocated storage space.
        %-----------------------------------------------------------------------        
        function [pq_new] = purge( pq )
            pq_new = PriorityQueue( pq.capacity );
            i = pq.top;
            while isfinite(i)
                pq_new = insert( pq_new, pq.eventArray{i} );
                i = pq.next(i);
            end
        end

        %-----------------------------------------------------------------------
        % Expand the capacity of the priority queue.
        %-----------------------------------------------------------------------        
        function [pq_new] = expand( pq )
            pq_new = PriorityQueue( 2*pq.capacity + 1 );
            
            pq_new.top                   = pq.top;
            pq_new.count                 = pq.count;
            pq_new.next(1:pq.capacity)   = pq.next;
            
            for j = 1:pq.capacity
                pq_new.eventArray{j} = pq.eventArray{j};
            end
        end
        
        %-----------------------------------------------------------------------
        % Dump the contents of the priority queue out to the previously opened 
        % file identifed by the file id <fid>.  This dump includes both active 
        % and zombie events.
        %-----------------------------------------------------------------------        
        function dump( pq, fid )
            fprintf(fid,'\n---------------------\n');
            fprintf(fid, '%10s %d %d %d \n', 'QUEUE', pq.top, pq.count, pq.capacity);
            
            for j = 1:pq.count
                fprintf('%10d %5d', j, pq.next(j));
                dump( pq.eventArray{j}, fid );
            end
            fprintf(fid,'---------------------\n');                        
        end
        
        %-----------------------------------------------------------------------
        % Dump the active events.
        %-----------------------------------------------------------------------        
        function dumpactive( pq, fid )
            j = pq.top;
            while isfinite(j)
                dump( pq.eventArray{j}, fid );
                j = pq.next(j);
            end
        end
        
    end
end