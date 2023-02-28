abstract type AbstractMigration end

struct Migrate <: AbstractMigration
    kind::Symbol
end
