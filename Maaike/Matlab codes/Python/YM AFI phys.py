#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service", token = "YOUR-API-KEY")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "secondaryIdentifier", "symbol",
    "interactions.participant2.secondaryIdentifier",
    "interactions.participant2.symbol", "interactions.details.relationshipType"
)

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "Gene list for S. cerevisiae 9 Nov 2020 8.56", code="A")
query.add_constraint("interactions.details.relationshipType", "=", "physical", code="B")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")

for row in query.rows():
    print(row["secondaryIdentifier"], row["symbol"], \
        row["interactions.participant2.secondaryIdentifier"], \
        row["interactions.participant2.symbol"], row["interactions.details.relationshipType"])
