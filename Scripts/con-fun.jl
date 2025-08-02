function convert_blanks_to_missing(df::DataFrame)
    for col in names(df)
        if eltype(df[!, col]) <: AbstractString
            df[!, col] = ifelse.(strip.(df[!, col]) .== "", missing, df[!, col])
        end
    end
    return df
end

function parse_datetime_column(dtc::AbstractString; format::Union{String, Nothing}=nothing)
    dtc_str = strip(dtc)
    
    if !isnothing(format)
        try
            return DateTime(dtc_str, format)
        catch e
            throw(ArgumentError("Failed to parse '$dtc_str' with provided format '$format': $(e.msg)"))
        end
    end
    
    formats = [
        "yyyy-mm-ddTHH:MM:SS",
        "yyyy-mm-dd HH:MM:SS", 
        "yyyy-mm-dd",
        "dd/mm/yyyy HH:MM:SS",
        "dd/mm/yyyy",
        "mm/dd/yyyy HH:MM:SS",
        "mm/dd/yyyy"
    ]
    
    for fmt in formats
        try
            return DateTime(dtc_str, fmt)
        catch
            continue
        end
    end
    
    throw(ArgumentError("Unable to parse datetime string '$dtc_str'. Tried formats: $(join(formats, ", ")). Consider providing a specific format using the 'format' keyword argument."))
end





function expand_doses_rows(ex_row::DataFrameRow, hours_between_doses)
    # 1) guard against any missing or invalid inputs
    if ismissing(hours_between_doses) ||
       ismissing(ex_row.ASTDTM)     ||
       ismissing(ex_row.AENDTM)     ||
       hours_between_doses ≤ 0      ||
       ex_row.AENDTM < ex_row.ASTDTM

        # just emit the row unchanged (as a NamedTuple) so nothing is lost
        return [NamedTuple(ex_row)]
    end

    # 2) compute how many doses
    duration_ms = Dates.value(ex_row.AENDTM - ex_row.ASTDTM)
    interval_ms = (86400000 * hours_between_doses) / 24
    n_doses    = Int(floor(duration_ms / interval_ms)) + 1

    # 3) build a NamedTuple for each dose‑time
    row_nt = NamedTuple(ex_row)
    out = Vector{NamedTuple}(undef, n_doses)
    for i in 1:n_doses
        # compute the fields that change
        newADTM    = ex_row.ASTDTM + Hour(hours_between_doses * (i - 1))
        newNFRLT   = ex_row.NFRLT     + hours_between_doses * (i - 1)
        newAVISITN = (newNFRLT ÷ 24) + 1

        # merge with the rest of the row
        merged = merge(row_nt, (
            ADTM       = newADTM,
            NFRLT      = newNFRLT,
            AVISITN    = newAVISITN,
            DOSE_INDEX = i
        ))

        out[i] = merged
    end

    return out
end



function expand_dose_dataframe(ex_df::DataFrame)
    # collect *all* expanded NamedTuples
    all_nts = NamedTuple[]
    for r in eachrow(ex_df)
        append!(all_nts, expand_doses_rows(r, r.HOURS_BETWEEN_DOSES))
    end
    # build one DataFrame — columns will be exactly those in ex_df plus DOSE_INDEX
    return DataFrame(all_nts)
end

function find_previous_dose(adpc::DataFrame, ex_exp::DataFrame)
    adpc.USUBJID = String.(adpc.USUBJID)
    ex_exp.USUBJID = String.(ex_exp.USUBJID)

    # Group ex_exp by USUBJID
    ex_by_usubjid = groupby(ex_exp, :USUBJID)

    # Prepare result holder
    prev_rows = Vector{Dict{Symbol, Any}}(undef, nrow(adpc))

    # Iterate over adpc
    for (i, row) in enumerate(eachrow(adpc))
        usubjid = row.USUBJID
        curr_adtm = row.ADTM
        group_key = (usubjid,)

        if haskey(ex_by_usubjid, group_key)
            subject_doses = ex_by_usubjid[group_key]
            # Filter rows where ADTM < current ADTM
            previous_doses = subject_doses[subject_doses.ADTM .< curr_adtm, :]

            if !isempty(previous_doses)
                # Take the row with latest ADTM
                idx = argmax(previous_doses.ADTM)
                last_dose = previous_doses[idx, :]

                prev_rows[i] = Dict(
                    :ADTM_prev => last_dose[:ADTM],
                    :EXDOSE_prev => last_dose[:EXDOSE],
                    :AVISIT_prev => last_dose[:AVISIT],
                    :AENDTM_prev => last_dose[:AENDTM]
                )
            else
                prev_rows[i] = Dict(
                    :ADTM_prev => missing,
                    :EXDOSE_prev => missing,
                    :AVISIT_prev => missing,
                    :AENDTM_prev => missing
                )
            end
        else
            prev_rows[i] = Dict(
                :ADTM_prev => missing,
                :EXDOSE_prev => missing,
                :AVISIT_prev => missing,
                :AENDTM_prev => missing
            )
        end
    end

    return hcat(adpc, DataFrame(prev_rows))
end

function find_next_dose(adpc::DataFrame, ex_exp::DataFrame)
    adpc.USUBJID = String.(adpc.USUBJID)
    ex_exp.USUBJID = String.(ex_exp.USUBJID)

    # Group ex_exp by USUBJID
    ex_by_usubjid = groupby(ex_exp, :USUBJID)

    # Prepare result holder
    next_rows = Vector{Dict{Symbol, Any}}(undef, nrow(adpc))

    # Iterate over adpc
    for (i, row) in enumerate(eachrow(adpc))
        usubjid = row.USUBJID
        curr_adtm = row.ADTM
        group_key = (usubjid,)

        if haskey(ex_by_usubjid, group_key)
            subject_doses = ex_by_usubjid[group_key]
            # Filter rows where ADTM >= current ADTM (next or concurrent dose)
            future_doses = subject_doses[subject_doses.ADTM .>= curr_adtm, :]

            if !isempty(future_doses)
                # Take the row with earliest ADTM (i.e., next dose)
                idx = argmin(future_doses.ADTM)
                next_dose = future_doses[idx, :]

                next_rows[i] = Dict(
                    :ADTM_next => next_dose[:ADTM],
                    :EXDOSE_next => next_dose[:EXDOSE],
                    :AVISIT_next => next_dose[:AVISIT],
                    :AENDTM_next => next_dose[:AENDTM]
                )
            else
                next_rows[i] = Dict(
                    :ADTM_next => missing,
                    :EXDOSE_next => missing,
                    :AVISIT_next => missing,
                    :AENDTM_next => missing
                )
            end
        else
            next_rows[i] = Dict(
                :ADTM_next => missing,
                :EXDOSE_next => missing,
                :AVISIT_next => missing,
                :AENDTM_next => missing
            )
        end
    end

    return hcat(adpc, DataFrame(next_rows))
end


function find_previous_nominal_time(adpc::DataFrame, ex_exp::DataFrame)
    adpc.USUBJID = String.(adpc.USUBJID)
    ex_exp.USUBJID = String.(ex_exp.USUBJID)

    ex_by_usubjid = groupby(ex_exp, :USUBJID)
    prev_nom_rows = Vector{Dict{Symbol, Any}}(undef, nrow(adpc))

    for (i, row) in enumerate(eachrow(adpc))
        usubjid = row.USUBJID
        nfrlt = row.NFRLT
        group_key = (usubjid,)

        if haskey(ex_by_usubjid, group_key) && !ismissing(nfrlt)
            subject_doses = ex_by_usubjid[group_key]

            # Only keep those with NFRLT < current
            previous_nominal = subject_doses[.!ismissing.(subject_doses.NFRLT) .&& subject_doses.NFRLT .< nfrlt, :]

            if !isempty(previous_nominal)
                idx = argmax(previous_nominal.NFRLT)
                last_nom = previous_nominal[idx, :]
                prev_nom_rows[i] = Dict(:NFRLT_prev => last_nom[:NFRLT])
            else
                prev_nom_rows[i] = Dict(:NFRLT_prev => missing)
            end
        else
            prev_nom_rows[i] = Dict(:NFRLT_prev => missing)
        end
    end

    return hcat(adpc, DataFrame(prev_nom_rows))
end


function find_next_nominal_time(adpc::DataFrame, ex_exp::DataFrame)
    adpc.USUBJID = String.(adpc.USUBJID)
    ex_exp.USUBJID = String.(ex_exp.USUBJID)

    ex_by_usubjid = groupby(ex_exp, :USUBJID)
    next_nom_rows = Vector{Dict{Symbol, Any}}(undef, nrow(adpc))

    for (i, row) in enumerate(eachrow(adpc))
        usubjid = row.USUBJID
        nfrlt = row.NFRLT
        group_key = (usubjid,)

        if haskey(ex_by_usubjid, group_key) && !ismissing(nfrlt)
            subject_doses = ex_by_usubjid[group_key]

            # Only keep those with NFRLT >= current
            future_nominal = subject_doses[.!ismissing.(subject_doses.NFRLT) .&& subject_doses.NFRLT .>= nfrlt, :]

            if !isempty(future_nominal)
                idx = argmin(future_nominal.NFRLT)
                next_nom = future_nominal[idx, :]
                next_nom_rows[i] = Dict(:NFRLT_next => next_nom[:NFRLT])
            else
                next_nom_rows[i] = Dict(:NFRLT_next => missing)
            end
        else
            next_nom_rows[i] = Dict(:NFRLT_next => missing)
        end
    end

    return hcat(adpc, DataFrame(next_nom_rows))
end