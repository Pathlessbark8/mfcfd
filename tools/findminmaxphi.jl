using ProgressMeter

function main()
    file_name = string(ARGS[1])
    data1 = read(file_name, String)
    splitdata = split(data1, "\n")
    numPoints = length(splitdata) - 1
    println("numPoints is ", numPoints)
    phi1_min, phi1_min_store, phi1_max, phi1_max_store = 1,1,1,0
    phi2_min, phi2_min_store, phi2_max, phi2_max_store = 1,1,1,0
    splitdata = @view split(data1, "\n")[2:end-1]
    @showprogress 1 "Computing ReadFile" for (idx, itm) in enumerate(splitdata)
        itmdata = split(itm)

        for vector_iter in 1:4
            store = parse(Float64,itmdata[vector_iter])
            if store > phi1_max_store
                phi1_max_store = store
                phi1_max = idx
            elseif store < phi1_min_store
                phi1_min_store = store
                phi1_min = idx
            end
        end
        
        for vector_iter in 5:8
            store = parse(Float64,itmdata[vector_iter])
            if store > phi2_max_store
                phi2_max_store = store
                phi2_max = idx
            elseif store < phi2_min_store
                phi2_min_store = store
                phi2_min = idx
            end
        end
    end
    println("Max Φ1 > ", phi1_max_store, " | Point Number - ", phi1_max)
    println("Min Φ1 > ", phi1_min_store, " | Point Number - ", phi1_min)
    println("Max Φ2 > ", phi2_max_store, " | Point Number - ", phi2_max)
    println("Min Φ2 > ", phi2_min_store, " | Point Number - ", phi2_min)
end

main()