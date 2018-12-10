%% add carnotUIBK (DO NOT MODIFY)
addpath('../../carnotUIBK_v1.3.2')

clear
clc

%% load the building with a name (MODIFY)
name_build_obj = 'building_EXAMPLE_PHPP_test_06_12_2018_v2';
[PLOT NUMBEROFZONES NUMBEROFINTERSECTIONS NUMBEROFWALLSINZONES NUMBEROFWINDOWSINZONES NUMBEROFWALLSININTERSECTIONS NUMBEROFWINDOWSINZONESANDINTERSECTIONS NUMBEROFGAINSINZONES NUMBEROFCONTROLS building] = load_building_automatic(name_build_obj);

%% parameters which allows you to automatize a run (MODIFY)
MODIFY = 1;
SIMULATE = 0;
ANALYSIS = 0;

%% run4simulink (DO NOT MODIFY)
run4simulink

%% MODIFY AND CREATE a building
if MODIFY
    % name of the EXCEL and XML file which the informations comes from (MODIFY)
    name_EXCEL = 'Copy_5_of_building_TEMPLATE.xlsx';
    name_XML = 'EXAMPLE.xml';
    name_PHPP = 'PHPP_EXAMPLE_rhocp.xlsx';
    
    % choose variant numbers (MODIFY)
    % modify:   -1     ... for no change
    %           0      ... for new variant with automatic number
    %           number ... for variant update or overwrite with new variant
    modify.geometry = 0;
    modify.construction = 0;
    modify.boundary = 0;
    modify.gains = 0;
    modify.hvac = 0;
    modify.thermalzone = 0;
    
    % choose the import mode (MODIFY)
    % import_mode: 'gbXML' ... import of gbXML file
    %              'excel' ... import of Excel file
    %              'PHPP'  ... import of PHPP !ATTENTION!: write also the name of the excel (it should be the excel template) you want that the data will be written on!
    import_mode = 'PHPP';
    
    % creates or modies the variants (DO NOT MODIFY)
    building = create_building_automatic(building, modify, import_mode, name_XML, name_EXCEL, name_PHPP, NUMBEROFZONES);
    run4simulink
end

%% SIMULATE a building
if SIMULATE
    % choose the variants to simulate (MODIFY)
    building.variant_geometry = 1;      % Excel
    building.variant_construction = 1;  % Excel
    building.variant_gains = 1;         % Excel
    building.variant_boundary = 1;      % Excel
    building.variant_thermalzone = 1;   % Excel
    building.variant_hvac = 1;          % Excel: with PKOM4 und Kaskadenlüftung
    
    % change model specific parameters if wished (MODIFY)
    building.preruntime = 30*24*3600;
    building.model = 'building_EXAMPLE_v1.slx';
    
    % load parameter for HVAC (MODIFY)
    HVAC_EXAMPLE = load('HVAC_EXAMPLE.mat');
    building.hvac(1) = HVAC_EXAMPLE.HVAC_EXAMPLE;
    
    % set some parameters of the simulations if wished (MODIFY)
    
%     theta_room = 21.0;
%     theta_room_neighbour = 21.0 - 1.0;
%     building.hvac(1).system(3).parameter.zone_1.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_2.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_3.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_4.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_5.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_6.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_7.teta_r_t = theta_room;
%     building.hvac(1).system(3).parameter.zone_8.teta_r_t = theta_room_neighbour;
%     building.hvac(1).system(3).parameter.zone_9.teta_r_t = theta_room_neighbour;
    
    % always use run4simulink before opening or running a model (DO NOT MODIFY)
    run4simulink
    
    % open the simulink model
    % building_EXAMPLE_v1
    
    % run a simulation (run directly and save automatically; or manually save it afterwards)
    building = building.simulate(1, 'description', 0);
    building = building.add_simulation(1, 'Innsbruck - ERV, 120m3/h, 230 g/h', saveAIB, saveBDB, saveBOUNDARY, saveHVAC, 1);
end

%% ANALYSIS of a simulation run
if ANALYSIS
    
end
