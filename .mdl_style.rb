all
exclude_rule 'MD033' # Allow inline HTML
exclude_rule 'MD041' # Do not require first line to be top level header
rule 'MD013', :line_length => 120 # Set line max length to 120
rule 'MD013', :tables => false # Disable line length on tables
