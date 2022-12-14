function [data] = interpolate_nan(data)

    if isnan(data(1))
        i = 2;

        while not isnan(data(i))
            i = i + 1;
        end

        data(1) = data(i);
    end

    if isnan(data(end))
        i = length(data);

        while ~isnan(data(i))
            i = i - 1;
        end

        data(1) = data(i);
    end

    for i = 2:length(data) - 1

        if isnan(data(i))

            if ~isnan(data(i + 1))
                data(i) = (data(i - 1) + data(i + 1)) / 2;
            else
                j = i + 2;

                while isnan(data(j))
                    j = j + 1;
                end

                data(i) = (data(i - 1) + data(i + j)) / 2;
            end

        end

    end

end
