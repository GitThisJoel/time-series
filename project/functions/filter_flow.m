function [flow] = fliter_flow(flow)
    m = mean(flow);
    [ci1, ~, ~] = ci(flow, 0.99);
    errs = flow < ci1;

    if errs(1)
        i = 2;

        while errs(i)
            i = i + 1;
        end

        flow(1) = flow(i);
    elseif errs(end)
        i = length(flow) - 1;

        while errs(i)
            i = i - 1;
        end

        flow(end) = flow(i);
    end

    for i = 2:length(errs) - 1

        if errs(i)

            if ~errs(i + 1)
                flow(i) = (flow(i - 1) + flow(i + 1)) / 2;
            else
                j = i + 2;

                while errs(j)
                    j = j + 1;
                end

                flow(i) = (flow(i - 1) + flow(j)) / 2;
            end

        end

    end

end
