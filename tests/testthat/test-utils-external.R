test_that("check_suggested_pkg errors when package is missing", {
  expect_error(
    check_suggested_pkg("definitely_not_a_real_pkg_xyz"),
    "definitely_not_a_real_pkg_xyz"
  )
})

test_that("check_suggested_pkg succeeds for base package", {
  expect_true(check_suggested_pkg("utils"))
})

test_that("check_external_tool errors when missing", {
  expect_error(
    check_external_tool("definitely_not_on_path_xyz"),
    "definitely_not_on_path_xyz"
  )
})

test_that("check_external_tool finds a known executable", {
  tool <- if (.Platform$OS.type == "windows") "cmd.exe" else "ls"
  path <- check_external_tool(tool)
  expect_true(nzchar(path))
})
