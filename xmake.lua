-- add rules: debug/release
add_rules("mode.debug", "mode.release")

-- define target
target("uwpflow")

    -- set kind
    set_kind("binary")

    -- add files
    add_headerfiles("Source/*.h")
    add_includedirs("Source")
    add_files("Source/*.cpp")

    -- set run envs
    set_rundir("Examples")
