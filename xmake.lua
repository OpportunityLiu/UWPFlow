-- add rules: debug/release
add_rules("mode.debug", "mode.release")

-- define target
target("uwpflow")

    -- set kind
    set_kind("binary")

    set_languages("cxx17")

    add_defines("WINDOWS", "ANSIPROTO")

    -- add files
    add_headerfiles("Source/*.h")
    add_includedirs("Source")
    add_files("Source/*.cpp")

    -- set run envs
    set_rundir("Examples")
