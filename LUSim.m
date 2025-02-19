classdef LUSim
% This class contains all the necessary pieces for simulating lattice
% ultrasensitivity
    properties
        param %struct containing relevant simulation parameters
        results %struct containing output data of simulations   
    end
    methods
        function obj = LUSim(param,results)
            obj.param = param;
            if nargin == 2 % makes it easy to initialize a new sim from the endstate of a previous one
                obj.param.M0 = results.Mend;
                obj.param.A0 = results.Aend;
                obj.param.cExp0 = results.cExp(end);
            end
            
        end
        
        function obj = runSim(obj) % runs simulation
            if obj.param.forScaling == 1
                obj.results = LUSimulatorForScaling(obj.param);
            else
                obj.results = LUSimulator(obj.param);
            end
        end
                
        function [] = plotActivity(obj) % plots simulated activity
            figure()
            plot(obj.results.ts/obj.param.nu,obj.results.As)
            xlabel('Time (s)')
            ylabel('CheA Activity')
        end
    end
end