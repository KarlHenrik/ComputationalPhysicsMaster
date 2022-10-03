struct MiniSort
    n::Int64
    x_sort::Vector{Float64}
    p::Vector{Int64}
    p_rev::Vector{Int64}

    x_sort_try::Vector{Float64}
    function MiniSort(positions)
        n = length(positions)
        p = sortperm(positions)
        x_sort = positions[p]
        p_rev = sortperm(p)
        x_sort_try = zero(positions)

        return new(n, x_sort, p, p_rev, x_sort_try)
    end
end

# Returns the sorted x with the new position, and the index of the sorted array
function try_sort!(minisort, new_idx, new_pos)
    (; n, x_sort, x_sort_try, p) = minisort
    i = p[new_idx]
    old_pos = x_sort[i]
    x_sort_try .= x_sort
    x_sort_try[i] = new_pos

    if new_pos > old_pos
        if i != n
            i += 1
            @inbounds while new_pos > x_sort_try[i]
                x_sort_try[i], x_sort_try[i - 1] = x_sort_try[i - 1], x_sort_try[i]
                if i == n
                    break
                end
                i += 1
            end
        end
    elseif new_pos < old_pos
        if i != 1
            i -= 1
            @inbounds while new_pos < x_sort_try[i]
                x_sort_try[i], x_sort_try[i + 1] = x_sort_try[i + 1], x_sort_try[i]
                if i == 1
                    break
                end 
                i -= 1
            end
        end
    end
   
    return x_sort_try, i
end

# Updates the positions, and therefore
function update_sort!(minisort, new_idx, new_pos)
    (; n, x_sort, p, p_rev) = minisort
    
    @inbounds sort_idx = p_rev[new_idx]
    @inbounds old_pos = x_sort[sort_idx]
    @inbounds x_sort[sort_idx] = new_pos

    if new_pos > old_pos
        if sort_idx != n
            r = sort_idx + 1
            @inbounds while new_pos > x_sort[r]
                x_sort[r], x_sort[r - 1] = x_sort[r - 1], x_sort[r]
                p[r], p[r - 1] = p[r - 1], p[r]
                p_rev[p[r]], p_rev[p[r - 1]] = p_rev[p[r - 1]] , p_rev[p[r]]

                r += 1
                if r == n + 1
                    break
                end
            end
        end
    elseif new_pos < old_pos
        if sort_idx != 1
            l = sort_idx - 1
            @inbounds while new_pos < x_sort[l]
                x_sort[l], x_sort[l + 1] = x_sort[l + 1], x_sort[l]
                p[l], p[l + 1] = p[l + 1], p[l]
                p_rev[p[l]], p_rev[p[l + 1]] = p_rev[p[l + 1]] , p_rev[p[l]]

                l -= 1
                if l == 0
                    break
                end 
            end
        end
    end
    return minisort
end