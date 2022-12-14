function [halt_ing_rep] = create_input_signal(raw_data)
    % FUNCTION create_input creates the input signal for part B

    % extract signals in raw data
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

    % construct signals of interest (higher abstraction level)
    flow_ing = flow_ing_1 + flow_ing_2; % [ton/h], kanske konvertera till m3/s?
    tol = 0.5;
    m = min(flow_ing);

    for i = 2:length(flow_ing) - 1

        if (abs(flow_ing(i)) < tol) || isnan(flow_ing(i))
            flow_ing(i) = (flow_ing(i - 1) + flow_ing(i + 1)) / 2;
        end

    end

%     disp("nan in flow ing " + sum(isnan(flow_ing)))
%     disp("ind ")

    % halt_ing_flot = (halt_ing_flot_1 .* flow_ing_1 + halt_ing_flot_2 .* flow_ing_2) ./ flow_ing;
    halt_raw_konc = (halt_raw_konc_1 .* flow_ing_1 + halt_raw_konc_2 .* flow_ing_2) ./ flow_ing;
    flow_raw_remill = flow_raw_remill_1 + flow_raw_remill_2;
    halt_scav_konc = (halt_scav_konc_1 .* flow_ing_1 + halt_scav_konc_2 .* flow_ing_2) ./ flow_ing;
    flow_scav_remill = flow_scav_remill_1 + flow_scav_remill_2;
    flow_ing_rep = (flow_raw_remill + flow_scav_remill);

    for i = 2:length(flow_ing_rep) - 1

        if (abs(flow_ing_rep(i)) < tol) || isnan(flow_ing_rep(i))
            flow_ing_rep(i) = (flow_ing_rep(i - 1) + flow_ing_rep(i + 1)) / 2;
        end

    end

    halt_ing_rep = (halt_raw_konc .* flow_raw_remill + halt_scav_konc .* flow_scav_remill) ./ flow_ing_rep;

    for i = 2:length(halt_ing_rep) - 1

        if (abs(halt_ing_rep(i)) < tol) || isnan(halt_ing_rep(i))
            halt_ing_rep(i) = (halt_ing_rep(i - 1) + halt_ing_rep(i + 1)) / 2;
        end

    end
    
    a = sum(isnan(halt_raw_konc)) % many nans
    b = sum(isnan(flow_raw_remill))
    c = sum(isnan(halt_scav_konc)) % many nans
    d = sum(isnan(flow_scav_remill))
    e = sum(isnan(flow_ing_rep))
    f = sum(isnan(halt_ing_rep))
    [a b c d e f]
end
