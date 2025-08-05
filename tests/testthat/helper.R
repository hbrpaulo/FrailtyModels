pkg_root <- normalizePath(file.path('..', '..'), mustWork = TRUE)
for (f in list.files(file.path(pkg_root, 'R'), pattern = '\\.[rR]$', full.names = TRUE)) {
  source(f)
}
