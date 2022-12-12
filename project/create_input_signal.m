function [halt_ing_rep] = create_input_signal()
    % FUNCTION create_input creates the input signal for part B

    load('raw_data.mat')

    % halt_konc = raw_data(:, 1);
    % halt_ing_flot_1 = raw_data(:, 2);
    % halt_ing_flot_2 = raw_data(:, 3);
    halt_raw_konc_1 = raw_data(:, 4);
    halt_raw_konc_2 = raw_data(:, 5);
    halt_scav_konc_1 = raw_data(:, 6);
    halt_scav_konc_2 = raw_data(:, 7);
    flow_raw_remill_1 = raw_data(:, 8);
    flow_raw_remill_2 = raw_data(:, 9);
    flow_scav_remill_1 = raw_data(:, 10);
    flow_scav_remill_2 = raw_data(:, 11);
    flow_ing_1 = raw_data(:, 13);
    flow_ing_2 = raw_data(:, 12);

    flow_ing = flow_ing_1 + flow_ing_2; % [ton/h], kanske konvertera till m3/s?
    % halt_ing_flot = (halt_ing_flot_1 .* flow_ing_1 + halt_ing_flot_2 .* flow_ing_2) ./ flow_ing;

    halt_raw_konc = (halt_raw_konc_1 .* flow_ing_1 + halt_raw_konc_2 .* flow_ing_2) ./ flow_ing;
    flow_raw_remill = flow_raw_remill_1 + flow_raw_remill_2;
    halt_scav_konc = (halt_scav_konc_1 .* flow_ing_1 + halt_scav_konc_2 .* flow_ing_2) ./ flow_ing;
    flow_scav_remill = flow_scav_remill_1 + flow_scav_remill_2;
    flow_ing_rep = (flow_raw_remill + flow_scav_remill);
    halt_ing_rep = (halt_raw_konc .* flow_raw_remill + halt_scav_konc .* flow_scav_remill) ./ flow_ing_rep; % %

end
