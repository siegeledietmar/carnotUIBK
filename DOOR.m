%% DOOR.m
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
classdef DOOR
    % DOOR
    
    properties
        name = '';
        X = 0;          % x coordinate left bottom
        Y = 0;          % y coordinate left bottom
        width = 0;
        height = 0;
        construction = '';
        model_heattrans = 0;
        model_cons = 0;
        view_factor = 1
        amb_factor = 1
        model_inf = 0;
        C = [0 0 0];
        n = [0 0 0];
        V = [0 0 0];
        control_i = 1;
        control_s = 1;
    end
    
    methods
        function obj = DOOR(name, X, Y, width, height, construction, model_heattrans, model_cons, view_factor, amb_factor, C, n, V, control_i, model_inf)
            if (nargin == 0)
            else
                obj.name = name;
                obj.X = X;
                obj.Y = Y;
                obj.width = width;
                obj.height = height;
                obj.construction = construction;
                obj.model_heattrans = model_heattrans;
                obj.model_cons = model_cons;
                obj.view_factor = view_factor;
                obj.amb_factor = amb_factor;
                obj.C = C;
                obj.n = n;
                obj.V = V;
                obj.control_i = control_i;
                obj.model_inf = model_inf;
            end
        end
        
        function area = calc_area(obj)
            % to calculate the area of the door
            area = obj.width * obj.height;
        end
    end
end

