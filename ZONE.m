%% ZONE.m
% ***********************************************************************
% This file is part of the uibkCARNOT Blockset.
% 
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%   Dietmar Siegele     dietmar.siegele@uibk.ac.at
%   Eleonora Leonardi   eleonora.leonardi@uibk.ac.at
% Additional Copyright for this file see list auf authors.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its 
%    contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
% THE POSSIBILITY OF SUCH DAMAGE.
% **********************************************************************
% 
%% carnotUIBK version 1.3
% Copyright (c) 2016-2018, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS,EL     2017-03-12   initial revision v1.0

%%
classdef ZONE
    % ZONE
    % heat flux always from 1 to 2!
    % 1 ... inside
    % 2 ... outside
    
    properties
        number = 0;
        name = '';                  % name of room
        rooms = {};                 % rooms
        model = -1;                 % 0 ... none/temperature, 1 ... ideal, 2 ... 1-node, 3 ... 2-node, 4 ... 2-node wo mass
        heated_area = 0;            % heated area [m²]
        heated_volume = 0;          % heated volume [m³]
        n50 = 0;                    % n50 value [1/h]
        t_ini = 20;                 % initial temperature [°C]
        phi_ini = 50;               % initial relative humidity [%]
        CO2_ini = 400;              % initial CO2 [ppm]
        VOC_ini = 0;                % inital VOC
        cp_spec = 204*3600;         % capacity of one-node [J/(K.m²)]
        timeconstant = 60;          % timeconstant [s]
        profile_timevalues = [0 1 8760]*3600; % time values for the profiles of Tc, Tr, rH, CO2, VOC, p; more than one setpoint has to be definied.. so here 1 hour of the year and then the rest of the year..
        profile_timeduration = [1 8760-1]*3600; % time duration for the profiles of Tc, Tr, rH, CO2, VOC, p (profiles if the model of the zone is 0 or 1)
        profile_Tc = [20 20];  % profile of Tc (profile to give if the model of the zone is 0 or 1)
        profile_Tr = [20 20];  % profile of Tr (profile to give if the model of the zone is 0 or 1)
        profile_rH = [50 50];  % profile of rH (profile to give if the model of the zone is 0 or 1)
        profile_CO2 = [400 400];  % profile of CO2 (profile to give if the model of the zone is 0 or 1)
        profile_VOC = [0 0];  % profile of VOC (profile to give if the model of the zone is 0 or 1)
        profile_p = [1e5 1e5];  % profile of p (profile to give if the model of the zone is 0 or 1)
        win_factor_seq = [];  % profile of fs during the year, time serie. One time serie for each window of the zone
        orientation_important = 1; % in the sum up of the walls for a room, if orientation_important = 1 the external orientation should be taken in consideration (walls with different orientations are seperated), if orientation_important = 0 the external orientation should not be taken in consideration (walls with different orientations are the same) 
       
        % matrix of the features of the walls and doors
        % 1...construction
        % 2...area
        % 3...boundary
        % 4...orientation_slope
        % 5...orientation_azimuth
        % 6...orientation_rotation
        % 7...model construction
        % 8...model heat transfer
        % 9...view factor
        % 10...ambient factor
        % 11...solar ratio
        % 12...number of boundary
        % 13...model_inf
        % 14...height
        % 15...C
        % 16...n
        % 17...V
        % 18...control_i
        matrix_wd = [];
        
        % matrix of the features of the windows
        % 1...construction
        % 2...width (sum of the widths when more windows)
        % 3...height (mean of the heights when more windows)
        % 4...boundary
        % 5...orientation_slope
        % 6...orientation_azimuth
        % 7...orientation_rotation
        % 8...model construction
        % 9...frame ratio
        % 10...glass area
        % 11...view factor
        % 12...ambient factor
        % 13...fs time (time serie)
        % 14...fs values (time serie)
        % 15...fd
        % 16...l*psi_installation
		% 17...length of glass
        % 18...control_s
        % 19...model_inf
        % 20...z (height of window above ground)
        % 21...C matrix [C0.0 C0.5 C1.0]
        % 22...n matrix [n0.0 n0.5 n1.0]
        % 23...V [door_closed door_open]
        % 24...control_i
        % 25...width (vector with the widths of the window when more window)
        % 26...height (vector with the heights of the window when more window)
        % 27...shadingtop (matrix with the two factors for each window)
        % 28...shadingleft (matrix with the two factors for each window)
        % 29...shadingright (matrix with the two factors for each window)
        % 30...shadinghorizont (matrix with the two factors for each window)
        matrix_wi = []; %
        
        % matrix of the features of the gains
        % 1...name
        % 2...model
        % 3...type
        % 4...timevalues
        % 5...timeduration
        % 6...values 1 (first serie)
        % 7...values 2 (first serie)
        % 8...values 3 (first serie)
        % 9...control
        matrix_gains = []; %
        
    end
    
    methods
        function obj = ZONE(number, name, rooms, model, orientation_important, t_ini, phi_ini, CO2_ini, VOC_ini, cp_spec, timeconstant, profile_Tc, profile_Tr, profile_time, profile_rH, profile_CO2, profile_VOC, profile_p)
            if (nargin == 4)
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
            elseif (nargin == 8)
                profile_Tc = t_ini;
                profile_Tr = phi_ini;
                profile_time = CO2_ini;
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
                obj.profile_Tc = profile_Tc;
                obj.profile_Tr = profile_Tr;
                obj.profile_timevalues = profile_time;
                obj.profile_timeduration = diff(profile_time);
            elseif (nargin == 12)
                profile_Tc = t_ini;
                profile_Tr = phi_ini;
                profile_time = CO2_ini;
                obj.profile_rH = VOC_ini;
                obj.profile_CO2 = cp_spec;
                obj.profile_VOC = timeconstant;
                obj.profile_p = profile_Tc;
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
                obj.profile_Tc = profile_Tc;
                obj.profile_Tr = profile_Tr;
                obj.profile_timevalues = profile_time;
                obj.profile_timeduration = diff(profile_time);
            elseif (nargin == 5)
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
            elseif (nargin == 11)
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
                obj.t_ini = t_ini;
                obj.phi_ini = phi_ini;
                obj.CO2_ini = CO2_ini;
                obj.VOC_ini = VOC_ini;
                obj.cp_spec = cp_spec;
                obj.timeconstant = timeconstant;
            elseif (nargin == 14)
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
                obj.t_ini = t_ini;
                obj.phi_ini = phi_ini;
                obj.CO2_ini = CO2_ini;
                obj.VOC_ini = VOC_ini;
                obj.cp_spec = cp_spec;
                obj.timeconstant = timeconstant;
                obj.profile_Tc = profile_Tc;
                obj.profile_Tr = profile_Tr;
                obj.profile_timevalues = profile_time;
                obj.profile_timeduration = diff(profile_time);
            elseif (nargin == 18)
                obj.number = number;
                obj.name = name;
                obj.rooms = rooms;
                obj.model = model;
                obj.orientation_important = orientation_important;
                obj.t_ini = t_ini;
                obj.phi_ini = phi_ini;
                obj.CO2_ini = CO2_ini;
                obj.VOC_ini = VOC_ini;
                obj.cp_spec = cp_spec;
                obj.timeconstant = timeconstant;
                obj.profile_Tc = profile_Tc;
                obj.profile_Tr = profile_Tr;
                obj.profile_timevalues = profile_time;
                obj.profile_timeduration = diff(profile_time);
                obj.profile_rH = profile_rH;
                obj.profile_CO2 = profile_CO2;
                obj.profile_VOC = profile_VOC;
                obj.profile_p = profile_p;
            else
                error('Not a valid ZONE, too less or too much parameters!')
            end
        end
        
    end
end