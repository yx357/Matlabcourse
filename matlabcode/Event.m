classdef Event
%-------------------------------------------------------------------------------
% Event
%
% Methods:
%   [e]    = Event( time, type, id1, id2 )
%   [flag] = isbefore( e1, e2 )
%   [flag] = involves( e, id )
%            dump( e, fid )
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
    properties
        time = nan;     % event time
        type = nan;     % event type
        id1  = nan;     % first disk id
        id2  = nan;     % second disk id
    end
    
    methods
        %-----------------------------------------------------------------------
        % Constructor.
        %-----------------------------------------------------------------------        
        function [e] = Event( time, type, id1, id2 )
            e.time = time;
            e.type = type;
            e.id1  = id1;
            e.id2  = id2;
        end
        
        %-----------------------------------------------------------------------
        % Returns true if event <e1> occurs before event <e2>.
        %-----------------------------------------------------------------------        
        function [flag] = isbefore( e1, e2 )
            flag = e1.time < e2.time;
        end

        %-----------------------------------------------------------------------
        % Returns true if event <e> involves the disk identified by <id>.
        %-----------------------------------------------------------------------        
        function [flag] = involves( e, id )
            flag = (e.id1 == id) | (e.id2 == id);
        end
        
        %-----------------------------------------------------------------------
        % Dump the event information out to the previously opened file 
        % identifed by the file id <fid>.
        %-----------------------------------------------------------------------
        function dump( e, fid )
            fprintf(fid,'%10s %10.2f %4d %4d %4d\n', 'EVENT', e.time, e.type, e.id1, e.id2 );
        end
    end
end