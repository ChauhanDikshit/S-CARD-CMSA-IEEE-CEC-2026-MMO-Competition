classdef Rotation<handle
    %% The class for string data for dependency (interactions among variable) of the problem"""
    properties (SetAccess=?ProblemMM)
        angle_max; % upper limit for the angles of rotations
        all_angles; % rotation angles for modes
        mat; % Rotation matrix for each basic function
    end
    methods 
        function output=Rotation() % requires functions ID and problem dimensionality 
            output.angle_max=pi; % upper limit for the angles of rotations
            output.all_angles=nan; % rotation angles for modes
            output.mat=nan; % Rotation matrix for each basic function
        end % function
    end % methods
end % class
