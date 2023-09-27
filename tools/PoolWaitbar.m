classdef PoolWaitbar < handle
    properties (SetAccess = immutable, GetAccess = private)
        Queue
        N
    end
    properties (Access = private, Transient)
        ClientHandle = []
        Count = 0
    end
    properties (SetAccess = immutable, GetAccess = private, Transient)
        Listener = []
    end
    methods (Access = private)
        function localIncrement(obj,message)
            obj.Count = 1 + obj.Count;
            waitbar(obj.Count/obj.N,obj.ClientHandle,... 
                {message; [num2str(round((obj.Count/obj.N)*100)) '$\%$ complete...']});
        end
    end
    methods
        function obj = PoolWaitbar(N, message)
            if nargin < 2
                message = 'PoolWaitbar';
            end
            obj.N = N;
            obj.ClientHandle = waitbar(0, 'Please wait...');
            obj.Queue = parallel.pool.DataQueue;
            obj.Listener = afterEach(obj.Queue, @(~) localIncrement(obj,message));
        end
        function increment(obj)
            send(obj.Queue, true);
        end
        function delete(obj)
            delete(obj.ClientHandle);
            delete(obj.Queue);
        end
    end
end

%% Developmental code
% 'Units','Normalized','Position', [0.5 0.5 0.25 0.1]