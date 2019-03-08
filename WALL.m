%% WALL.m
% ***********************************************************************
% This file is part of the uibkCARNOT Blockset.
% 
% Copyright (c) 2016-2019, University of Innsbruck, Unit for Energy 
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
%% carnotUIBK version 2.0
% Copyright (c) 2016-2019, University of Innsbruck, Unit for Energy 
% Efficient Building.
%
% Author    Date         Description
% DS,EL     2017-03-12   initial revision v1.0
% DS,EL	    2019-01-24   updates for GUI v2.0

%%
classdef WALL
    % WALL
    
    properties
        name = '';
        boundary = '';  % boundary object ROOM or AMBIENT, GROUND
        X = 0;          % x coordinate left bottom
        Y = 0;          % y coordinate left bottom
        Z = 0;          % z coordinate left bottom
        width = 0;
        height = 0;
        orientation_slope = 0;      % 0° ... ceiling, 90° ... vertical, 180° ... floor
        orientation_azimuth = 0;    % 0° ... south, -90° ... east
        orientation_rotation = 0;   % rotation
        inside = 1;                 % construction: 1 ... inside, 0 ... outside          
        construction = '';
        model_cons = 0;         % 0 ... UA, 1 ... RC, 2 ... hygrothermal
        model_heattrans = 0;    % 0 ... none, 1 ... const, 2 ... complex
        view_factor = 1;
        amb_factor = 1;
        windows = [];
        doors = [];
        model_inf = -1;
        C = [0 0 0];
        n = [0 0 0];
        V = [0 0 0];
        control_i = 1;
    end
    
    methods
        
        function obj = WALL(name, boundary, X, Y, Z, width, height, orientation_slope, orientation_azimuth, orientation_rotation, inside, construction, model_cons, model_heattrans, view_factor, amb_factor, model_inf, C, n, V, control_i, windows, doors)
            if (nargin == 0)
            else
                obj.name = name;
                obj.boundary = boundary;
                obj.X = X;
                obj.Y = Y;
                obj.Z = Z;
                obj.width = width;
                obj.height = height;
                obj.orientation_slope = orientation_slope;
                obj.orientation_azimuth = orientation_azimuth;
                obj.orientation_rotation = orientation_rotation;
                obj.inside = inside;
                obj.construction = construction;
                obj.model_cons = model_cons;
                obj.model_heattrans = model_heattrans;
                obj.view_factor = view_factor;
                obj.amb_factor = amb_factor;
                obj.windows = windows;
                obj.doors = doors; 
                obj.model_inf = model_inf;
                obj.C = C;
                obj.n = n;
                obj.V = V;
                obj.control_i = control_i;    
            end
            
        end
        
        function area = calc_area(obj)
            % to calcolate the area of a wall subtracting the area of the
            % windows and the doors
            area = obj.width * obj.height;
            for jj = 1:length(obj.windows)
                area = area - obj.windows(jj).calc_area();
            end
            for jj = 1:length(obj.doors)
                area = area - obj.doors(jj).calc_area();
            end
        end
        
        function plot(obj, color)
            % to plot the wall with window and door
            % 1 ... codex of the color of the wall
           
            if obj.X <= -999
                
            elseif length(obj.X) > 1
                x = obj.X;
                y = obj.Y;
                z = obj.Z;
                
                patch(x,y,z,color)
                
                x = obj.X(2);
                y = obj.Y(2);
                z = obj.Z(2);
                
                gam = obj.orientation_azimuth;
                al = obj.orientation_slope;
                be = obj.orientation_rotation;
                
                Rx = [1 0 0 0
                    0 cos(al*pi/180) -sin(al*pi/180) 0
                    0 sin(al*pi/180) cos(al*pi/180) 0
                    0 0 0 1];
                
                Ry = [cos(be*pi/180) 0 sin(be*pi/180) 0
                    0 1 0 0
                    -sin(be*pi/180) 0 cos(be*pi/180) 0
                    0 0 0 1];
                
                Rz = [cos(gam*pi/180) -sin(gam*pi/180) 0 0
                    sin(gam*pi/180) cos(gam*pi/180) 0 0
                    0 0 1 0
                    0 0 0 1];
                
                position_0 = [x
                    y
                    z
                    1];
                
                coo_1 = [];
                coo_2 = [];
                coo_3 = [];
                coo_4 = [];
                
                coo_1 = [0; 0; 0; 1];
                coo_2 = Rz*Rx*Ry*[obj.width; 0; 0; 1];
                coo_3 = Rz*Rx*Ry*[obj.width; obj.height; 0; 1];
                coo_4 = Rz*Rx*Ry*[0; obj.height; 0; 1];
                
                coo_1 = round(coo_1*1000)/1000 + position_0;
                coo_2 = round(coo_2*1000)/1000 + position_0;
                coo_3 = round(coo_3*1000)/1000 + position_0;
                coo_4 = round(coo_4*1000)/1000 + position_0;
                
                if ~isempty(obj.windows)
                    for ii=1:size(obj.windows,2)
                        XW = obj.windows(ii).X;
                        YW = obj.windows(ii).Y;
                        widthW = obj.windows(ii).width;
                        heightW = obj.windows(ii).height;

                        coo_1W = [];
                        coo_2W = [];
                        coo_3W = [];
                        coo_4W = [];

                        coo_1W = Rz*Rx*Ry*[XW; YW; 0; 1];
                        coo_2W = Rz*Rx*Ry*[XW+widthW; YW; 0; 1];
                        coo_3W = Rz*Rx*Ry*[XW+widthW; YW+heightW; 0; 1];
                        coo_4W = Rz*Rx*Ry*[XW; YW+heightW; 0; 1];

                        coo_1W = round(coo_1W*1000)/1000 + position_0;
                        coo_2W = round(coo_2W*1000)/1000 + position_0;
                        coo_3W = round(coo_3W*1000)/1000 + position_0;
                        coo_4W = round(coo_4W*1000)/1000 + position_0;

                        patch([coo_1W(1),coo_2W(1),coo_3W(1),coo_4W(1)],[coo_1W(2),coo_2W(2),coo_3W(2),coo_4W(2)],[coo_1W(3),coo_2W(3),coo_3W(3),coo_4W(3)],'cyan')
                    end
                end

                if ~isempty(obj.doors)
                    for ii=1:size(obj.doors,2)
                        XD = obj.doors(ii).X;
                        YD = obj.doors(ii).Y;
                        widthD = obj.doors(ii).width;
                        heightD = obj.doors(ii).height;

                        coo_1D = [];
                        coo_2D = [];
                        coo_3D = [];
                        coo_4D = [];

                        coo_1D = Rz*Rx*Ry*[XD; YD; 0; 1];
                        coo_2D = Rz*Rx*Ry*[XD+widthD; YD; 0; 1];
                        coo_3D = Rz*Rx*Ry*[XD+widthD; YD+heightD; 0; 1];
                        coo_4D = Rz*Rx*Ry*[XD; YD+heightD; 0; 1];

                        coo_1D = round(coo_1D*1000)/1000 + position_0;
                        coo_2D = round(coo_2D*1000)/1000 + position_0;
                        coo_3D = round(coo_3D*1000)/1000 + position_0;
                        coo_4D = round(coo_4D*1000)/1000 + position_0;

                        patch([coo_1D(1),coo_2D(1),coo_3D(1),coo_4D(1)],[coo_1D(2),coo_2D(2),coo_3D(2),coo_4D(2)],[coo_1D(3),coo_2D(3),coo_3D(3),coo_4D(3)], [0.7, 0.5, 0])
                    end
                end
            else
                x = obj.X;
                y = obj.Y;
                z = obj.Z;

                gam = obj.orientation_azimuth;
                al = obj.orientation_slope;
                be = obj.orientation_rotation;

                Rx = [1 0 0 0
                    0 cos(al*pi/180) -sin(al*pi/180) 0
                    0 sin(al*pi/180) cos(al*pi/180) 0
                    0 0 0 1];

                Ry = [cos(be*pi/180) 0 sin(be*pi/180) 0
                    0 1 0 0
                    -sin(be*pi/180) 0 cos(be*pi/180) 0
                    0 0 0 1];

                Rz = [cos(gam*pi/180) -sin(gam*pi/180) 0 0
                    sin(gam*pi/180) cos(gam*pi/180) 0 0
                    0 0 1 0
                    0 0 0 1];

                position_0 = [x
                    y
                    z
                    1];

                coo_1 = [];
                coo_2 = [];
                coo_3 = [];
                coo_4 = [];

                coo_1 = [0; 0; 0; 1];
                coo_2 = Rz*Rx*Ry*[obj.width; 0; 0; 1];
                coo_3 = Rz*Rx*Ry*[obj.width; obj.height; 0; 1];
                coo_4 = Rz*Rx*Ry*[0; obj.height; 0; 1];

                coo_1 = round(coo_1*1000)/1000 + position_0;
                coo_2 = round(coo_2*1000)/1000 + position_0;
                coo_3 = round(coo_3*1000)/1000 + position_0;
                coo_4 = round(coo_4*1000)/1000 + position_0;

                patch([coo_1(1),coo_2(1),coo_3(1),coo_4(1)],[coo_1(2),coo_2(2),coo_3(2),coo_4(2)],[coo_1(3),coo_2(3),coo_3(3),coo_4(3)],color)

                if ~isempty(obj.windows)
                    for ii=1:size(obj.windows,2)
                        XW = obj.windows(ii).X;
                        YW = obj.windows(ii).Y;
                        widthW = obj.windows(ii).width;
                        heightW = obj.windows(ii).height;

                        coo_1W = [];
                        coo_2W = [];
                        coo_3W = [];
                        coo_4W = [];

                        coo_1W = Rz*Rx*Ry*[XW; YW; 0; 1];
                        coo_2W = Rz*Rx*Ry*[XW+widthW; YW; 0; 1];
                        coo_3W = Rz*Rx*Ry*[XW+widthW; YW+heightW; 0; 1];
                        coo_4W = Rz*Rx*Ry*[XW; YW+heightW; 0; 1];

                        coo_1W = round(coo_1W*1000)/1000 + position_0;
                        coo_2W = round(coo_2W*1000)/1000 + position_0;
                        coo_3W = round(coo_3W*1000)/1000 + position_0;
                        coo_4W = round(coo_4W*1000)/1000 + position_0;

                        patch([coo_1W(1),coo_2W(1),coo_3W(1),coo_4W(1)],[coo_1W(2),coo_2W(2),coo_3W(2),coo_4W(2)],[coo_1W(3),coo_2W(3),coo_3W(3),coo_4W(3)],'cyan')
                    end
                end

                if ~isempty(obj.doors)
                    for ii=1:size(obj.doors,2)
                        XD = obj.doors(ii).X;
                        YD = obj.doors(ii).Y;
                        widthD = obj.doors(ii).width;
                        heightD = obj.doors(ii).height;

                        coo_1D = [];
                        coo_2D = [];
                        coo_3D = [];
                        coo_4D = [];

                        coo_1D = Rz*Rx*Ry*[XD; YD; 0; 1];
                        coo_2D = Rz*Rx*Ry*[XD+widthD; YD; 0; 1];
                        coo_3D = Rz*Rx*Ry*[XD+widthD; YD+heightD; 0; 1];
                        coo_4D = Rz*Rx*Ry*[XD; YD+heightD; 0; 1];

                        coo_1D = round(coo_1D*1000)/1000 + position_0;
                        coo_2D = round(coo_2D*1000)/1000 + position_0;
                        coo_3D = round(coo_3D*1000)/1000 + position_0;
                        coo_4D = round(coo_4D*1000)/1000 + position_0;

                        patch([coo_1D(1),coo_2D(1),coo_3D(1),coo_4D(1)],[coo_1D(2),coo_2D(2),coo_3D(2),coo_4D(2)],[coo_1D(3),coo_2D(3),coo_3D(3),coo_4D(3)], [0.7, 0.5, 0])
                    end
                end
            end
            view(3)
        end
    end
end

